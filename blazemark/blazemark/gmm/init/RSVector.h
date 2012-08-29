//=================================================================================================
/*!
//  \file blazemark/gmm/init/RSVector.h
//  \brief Header file for the GMM++ sparse vector initialization functions
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

#ifndef _BLAZEMARK_GMM_INIT_RSVECTOR_H_
#define _BLAZEMARK_GMM_INIT_RSVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <gmm/gmm_vector.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Indices.h>


namespace blazemark {

namespace gmm {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name GMM++ initialization functions */
//@{
template< typename Type >
void init( ::gmm::rsvector<Type>& v, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given sparse vector.
//
// \param v The sparse vector to be initialized.
// \param nonzeros The number of non-zero elements.
// \return void
//
// This function initializes the given sparse vector with random values. The vector will
// be filled with \a nonzeros non-zero elements, whose indices will be randomly determined.
*/
template< typename Type >  // Data type of the vector
void init( ::gmm::rsvector<Type>& v, size_t nonzeros )
{
   const size_t N( vect_size( v ) );

   ::blazemark::Indices indices( N, nonzeros );
   for( ::blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
      v[*it] = ::blaze::rand<Type>();
   }
}
//*************************************************************************************************

} // namespace gmm

} // namespace blazemark

#endif
