//=================================================================================================
/*!
//  \file blaze/util/policies/ConstantGrowth.h
//  \brief Header file for the ConstantGrowth policy classes.
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

#ifndef _BLAZE_UTIL_POLICIES_CONSTANTGROWTH_H_
#define _BLAZE_UTIL_POLICIES_CONSTANTGROWTH_H_


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
/*!\brief Constant growth policy class.
// \ingroup util
//
// The ConstantGrowth policy class implements a constant growth strategy. It can be customized
// for any purpose: the \a Growth template argument specifies the constant increase of the given
// size.
*/
template< size_t Growth >
struct ConstantGrowth
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
/*!\brief Specialization of the constant growth policy for 0 growth.
// \ingroup util
*/
template<>
struct ConstantGrowth<0>;
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
inline size_t ConstantGrowth<Growth>::operator()( size_t old, size_t minimum ) const
{
   const size_t needed( max<size_t>( old+Growth, minimum ) );
   return ( ( needed )?( 4 * ( (needed-1)/4+1 ) ):( 0 ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
