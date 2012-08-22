//=================================================================================================
/*!
//  \file blaze/util/policies/LinearGrowth.h
//  \brief Header file for the LinearGrowth policy classes.
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

#ifndef _BLAZE_UTIL_POLICIES_LINEARGROWTH_H_
#define _BLAZE_UTIL_POLICIES_LINEARGROWTH_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Functions.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Linear growth policy class.
// \ingroup util
//
// The LinearGrowth policy class implements a linear growth strategy. It can be customized for
// any purpose: the \a Growth template argument specifies the factor of the size growth.
*/
template< size_t Growth >
struct LinearGrowth
{
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t operator()( size_t oldSize, size_t minSize ) const;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the linear growth policy for 0 growth.
// \ingroup util
*/
template<>
struct LinearGrowth<0>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the linear growth policy for 1 growth.
// \ingroup util
*/
template<>
struct LinearGrowth<1>;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a new size depending on the given old size and the required minimum size.
//
// \param old The old size.
// \param minimum The required minimum size.
// \return The new size (at least the required minimum size).
*/
template< size_t Growth >
inline size_t LinearGrowth<Growth>::operator()( size_t old, size_t minimum ) const
{
   const size_t needed( max<size_t>( old*Growth, minimum ) );
   return ( ( needed )?( 4 * ( (needed-1)/4+1 ) ):( 0 ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
