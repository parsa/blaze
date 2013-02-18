//=================================================================================================
/*!
//  \file blaze/util/policies/Deallocate.h
//  \brief Header file for the Deallocate policy classes.
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

#ifndef _BLAZE_UTIL_POLICIES_DEALLOCATE_H_
#define _BLAZE_UTIL_POLICIES_DEALLOCATE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Memory.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Deallocate policy class.
// \ingroup util
//
// The Deallocate deletion policy is the according deletion policy for arrays allocated via
// the blaze::allocate function. It uses deallocate to free the resource. Note that the delete
// operation is NOT permitted for inclomplete types (i.e. declared but undefined data types).
// The attempt to apply a PtrDelete functor to a pointer to an object of incomplete type
// results in a compile time error!
*/
struct Deallocate
{
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename Type >
   inline void operator()( Type ptr ) const;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of the deallocate deletion policy.
//
// \param ptr The pointer to delete.
// \return void
//
// This function frees the given pointer resource via the blaze::deallocate function. Note that
// the delete operation is NOT permitted for incomplete types (i.e. declared but undefined data
// types). The attempt to use this function for a pointer to an object of incomplete type results
// in a compile time error!
*/
template< typename Type >
inline void Deallocate::operator()( Type ptr ) const
{
   deallocate( ptr );
}
//*************************************************************************************************

} // namespace blaze

#endif
