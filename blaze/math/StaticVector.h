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
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Rand();
   explicit inline Rand( Type min, Type max );
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator StaticVector<Type,N,TF>() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   StaticVector<Type,N,TF> vector_;  //!< The random vector.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default constructor of the Rand specialization for StaticVector.
*/
template< typename Type  // Data type of the vector
        , size_t N       // Number of elements
        , bool TF >      // Transpose flag
inline Rand< StaticVector<Type,N,TF> >::Rand()
   : vector_()  // The random vector
{
   for( size_t i=0UL; i<N; ++i ) {
      vector_[i] = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Range constructor of the Rand specialization for StaticVector.
//
// \param min The smallest possible value for a vector element.
// \param max The largest possible value for a vector element.
*/
template< typename Type  // Data type of the vector
        , size_t N       // Number of elements
        , bool TF >      // Transpose flag
inline Rand< StaticVector<Type,N,TF> >::Rand( Type min, Type max )
   : vector_()  // The random vector
{
   for( size_t i=0UL; i<N; ++i ) {
      vector_[i] = rand<Type>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion to the created random StaticVector.
//
// \return The random vector.
*/
template< typename Type  // Data type of the vector
        , size_t N       // Number of elements
        , bool TF >      // Transpose flag
inline Rand< StaticVector<Type,N,TF> >::operator StaticVector<Type,N,TF>() const
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
