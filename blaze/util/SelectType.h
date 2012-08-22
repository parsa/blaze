//=================================================================================================
/*!
//  \file blaze/util/SelectType.h
//  \brief Header file for the SelectType class template
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

#ifndef _BLAZE_UTIL_SELECTTYPE_H_
#define _BLAZE_UTIL_SELECTTYPE_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type selection.
// \ingroup util
//
// The SelectType class template selects one of the two given types \a T1 and \a T2 depending
// on the \a Select template argument. In case the \a Select compile time constant expression
// evaluates to \a true, the member type definition \a Type is set to \a T1. In case \a Select
// evaluates to \a false, \a Type is set to \a T2.
*/
template< bool Select    // Compile time selection
        , typename T1    // Type to be selected if Select=true
        , typename T2 >  // Type to be selected if Select=false
struct SelectType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef T1  Type;  //!< The selected type.
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the SwitchType class template.
// \ingroup util
//
// This specialization of the SelectType template is selected in case the \a Select compile time
// constant expression evaluates to \a false. The member type definition is set to the second
// given type \a T2.
*/
template< typename T1    // Type not to be selected
        , typename T2 >  // Type to be selected
struct SelectType<false,T1,T2>
{
 public:
   //**********************************************************************************************
   typedef T2  Type;  //!< The selected type.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
