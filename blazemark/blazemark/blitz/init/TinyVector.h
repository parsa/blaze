//=================================================================================================
/*!
//  \file blazemark/blitz/init/TinyVector.h
//  \brief Header file for the Blitz++ tiny vector initialization functions
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

#ifndef _BLAZEMARK_BLITZ_INIT_TINYVECTOR_H_
#define _BLAZEMARK_BLITZ_INIT_TINYVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blitz/tinyvec2.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace blitz {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Blitz++ initialization functions */
//@{
template< typename Type, int N >
void init( ::blitz::TinyVector<Type,N>& v );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given tiny vector.
//
// \param v The tiny vector to be initialized.
// \return void
//
// This function initializes the given tiny vector with random values.
*/
template< typename Type  // Data type of the vector
        , int N >        // Number of elements
void init( ::blitz::TinyVector<Type,N>& v )
{
   for( int i=0; i<N; ++i ) {
      v(i) = ::blaze::rand<Type>();
   }
}
//*************************************************************************************************

} // namespace blitz

} // namespace blazemark

#endif
