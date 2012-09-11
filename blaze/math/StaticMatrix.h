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

#ifndef _BLAZE_MATH_STATIC_MATRIX_H_
#define _BLAZE_MATH_STATIC_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/system/Precision.h>


namespace blaze {

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
