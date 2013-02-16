//=================================================================================================
/*!
//  \file blaze/util/policies/NoDelete.h
//  \brief Header file for the NoDelete policy classes.
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

#ifndef _BLAZE_UTIL_POLICIES_NODELETE_H_
#define _BLAZE_UTIL_POLICIES_NODELETE_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief No-delete policy class.
// \ingroup util
*/
template< typename Type >
struct NoDelete
{
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void operator()( const Type& ptr ) const;
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
/*!\brief Implementation of the no-delete policy.
//
// \param ptr The pointer to delete.
// \return void
*/
template< typename Type >
inline void NoDelete<Type>::operator()( const Type& /*ptr*/ ) const
{}
//*************************************************************************************************

} // namespace blaze

#endif
