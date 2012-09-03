//=================================================================================================
/*!
//  \file blazemark/flens/init/DenseVector.h
//  \brief Header file for the FLENS dense vector initialization functions
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

#ifndef _BLAZEMARK_FLENS_INIT_DENSEVECTOR_H_
#define _BLAZEMARK_FLENS_INIT_DENSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <flens/vectortypes/impl/densevector.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace flens {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name FLENS initialization functions */
//@{
template< typename Type >
void init( ::flens::DenseVector< ::flens::Array<Type> >& v );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given dense vector.
//
// \param v The dense vector to be initialized.
// \return void
//
// This function initializes the given dense vector with random values.
*/
template< typename Type >  // Data type of the vector
void init( ::flens::DenseVector< ::flens::Array<Type> >& v )
{
   typedef typename ::flens::DenseVector< ::flens::Array<Type> >::IndexType  IndexType;

   for( IndexType i=v.firstIndex(); i<=v.lastIndex(); ++i ) {
      v(i) = ::blaze::rand<Type>();
   }
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark

#endif
