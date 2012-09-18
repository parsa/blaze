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

#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/system/Precision.h>


namespace blaze {

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
