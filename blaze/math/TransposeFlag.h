//=================================================================================================
/*!
//  \file blaze/math/TransposeFlag.h
//  \brief Header file for the vector transpose flag types
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

#ifndef _BLAZE_MATH_TRANSPOSEFLAG_H_
#define _BLAZE_MATH_TRANSPOSEFLAG_H_


namespace blaze {

//=================================================================================================
//
//  VECTOR TRANSPOSE FLAG TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Transpose flag for column vectors.
//
// Via this flag it is possible to specify vectors as column vectors. The following example
// demonstrates the setup of a 3-dimensional column vector:

   \code
   using blaze::columnVector;
   blaze::StaticVector<int,3UL,columnVector> v( 1, 2, 3 );
   \endcode
*/
const bool columnVector = false;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transpose flag for row vectors.
//
// Via this flag it is possible to specify vectors as row vectors. The following example
// demonstrates the setup of a 3-dimensional row vector:

   \code
   using blaze::rowVector;
   blaze::StaticVector<int,3UL,rowVector> v( 1, 2, 3 );
   \endcode
*/
const bool rowVector = true;
//*************************************************************************************************

} // namespace blaze

#endif
