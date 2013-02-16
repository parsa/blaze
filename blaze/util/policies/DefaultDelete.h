//=================================================================================================
/*!
//  \file blaze/util/policies/DefaultDelete.h
//  \brief Header file for the DefaultDelete policy classes.
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

#ifndef _BLAZE_UTIL_POLICIES_DEFAULTDELETE_H_
#define _BLAZE_UTIL_POLICIES_DEFAULTDELETE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/checked_delete.hpp>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default C++ deletion policy class.
// \ingroup util
//
// The DefaultDelete deletion policy is the standard delete for resources allocated via the new
// operator. It uses delete or array delete (depending on the template argument) to free the
// resource:

   \code
   class Resource { ... };
   
   DefaultDelete<Resource> ptrDelete       // Uses delete to free resources
   DefaultDelete<Resource[]> arrayDelete;  // Uses array delete to free resources
   \endcode

// Note the explicit use of empty array bounds to configure DefaultDelete to use array delete
// instead of delete. Also note that the delete operation is NOT permitted for incomplete types
// (i.e. declared but undefined data types). The attempt to apply a DefaultDelete functor to a
// pointer or array to an object of incomplete type results in a compile time error!
*/
template< typename Type >
struct DefaultDelete
{
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void operator()( Type* ptr ) const;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Implementation of the default delete policy.
//
// \param ptr The pointer to delete.
// \return void
//
// This function frees the given pointer resource via delete. Note that the delete operation
// is NOT permitted for incomplete types (i.e. declared but undefined data types). The attempt
// to use this function for a pointer to an object of incomplete type results in a compile time
// error!
*/
template< typename Type >
inline void DefaultDelete<Type>::operator()( Type* ptr ) const
{
   boost::checked_delete( ptr );
}
//*************************************************************************************************




//=================================================================================================
//
//  SPECIALIZATION FOR ARRAYS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DefaultDelete class template for arrays.
// \ingroup util
//
// This specialization of the DefaultDelete class template uses array delete to free the
// allocated resource.
*/
template< typename Type >
struct DefaultDelete<Type[]>
{
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void operator()( Type* ptr ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the default array delete policy.
//
// \param ptr The pointer to delete.
// \return void
//
// This function frees the given array resource via array delete. Note that the delete operation
// is NOT permitted for incomplete types (i.e. declared but undefined data types). The attempt
// to use this function for a pointer to an object of incomplete type results in a compile time
// error!
*/
template< typename Type >
inline void DefaultDelete<Type[]>::operator()( Type* ptr ) const
{
   boost::checked_array_delete( ptr );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
