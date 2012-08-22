//=================================================================================================
/*!
//  \file blaze/util/mpl/IfNot.h
//  \brief Header file for the IfNot class template
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

#ifndef _BLAZE_UTIL_MPL_IFNOT_H_
#define _BLAZE_UTIL_MPL_IFNOT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/SelectType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type value comparison.
// \ingroup mpl
//
// The IfNot class template selects one of the two given types \a T2 and \a T3 depending on \a T1.
// In case \a T1::value evaluates to \a false, the member type definition \a Type is set to \a T2.
// In case \a T1::value evaluates to \a true, \a Type is set to \a T3.
*/
template< typename T1    // Type of the condition
        , typename T2    // Type to be selected if T1::value=false
        , typename T3 >  // Type to be selected if T1::value=true
struct IfNot
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< !T1::value, T2, T3 >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
