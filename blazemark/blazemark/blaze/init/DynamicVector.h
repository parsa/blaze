//=================================================================================================
/*!
//  \file blazemark/blaze/init/DynamicVector.h
//  \brief Header file for the Blaze dynamic vector initialization functions
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

#ifndef _BLAZEMARK_BLAZE_INIT_DYNAMICVECTOR_H_
#define _BLAZEMARK_BLAZE_INIT_DYNAMICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <vector>
#include <blaze/math/DynamicVector.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace blaze {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Blaze initialization functions */
//@{
template< typename Type, bool TF >
void init( ::blaze::DynamicVector<Type,TF>& v );

template< typename Type, bool TF >
void init( ::std::vector< ::blaze::DynamicVector<Type,TF> >& v );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given dynamic vector.
//
// \param v The dynamic vector to be initialized.
// \return void
//
// This function initializes the given dynamic vector with random values.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
void init( ::blaze::DynamicVector<Type,TF>& v )
{
   const size_t N( v.size() );

   for( size_t i=0UL; i<N; ++i ) {
      v[i] = ::blaze::rand<Type>( 0, 10 );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given vector of dynamic vector.
//
// \param v The vector of dynamic vectors to be initialized.
// \return void
//
// This function initializes all dynamic vectors in the given vector with random values.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
void init( ::std::vector< ::blaze::DynamicVector<Type,TF> >& v )
{
   const size_t size( v.size() );

   for( size_t i=0UL; i<size; ++i ) {
      init( v[i] );
   }
}
//*************************************************************************************************

} // namespace blaze

} // namespace blazemark

#endif
