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

#ifndef _BLAZE_MATH_STATIC_VECTOR_H_
#define _BLAZE_MATH_STATIC_VECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/dense/StaticVector.h>
#include <blaze/system/Precision.h>


namespace blaze {

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
