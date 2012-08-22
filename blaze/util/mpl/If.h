//=================================================================================================
/*!
//  \file blaze/util/mpl/If.h
//  \brief Header file for the If class template
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

#ifndef _BLAZE_UTIL_MPL_IF_H_
#define _BLAZE_UTIL_MPL_IF_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type selection.
// \ingroup mpl
//
// The IfTrue class template selects one of the two given types \a T1 and \a T2 depending
// on the \a Condition template argument. In case the \a Condition compile time constant
// expression evaluates to \a true, the member type definition \a Type is set to \a T1.
// In case \a Condition evaluates to \a false, \a Type is set to \a T2.
*/
template< bool Condition  // Compile time selection
        , typename T1     // Type to be selected if Condition=true
        , typename T2 >   // Type to be selected if Condition=false
struct IfTrue
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
/*!\brief Specialization of the IfTrue class template.
// \ingroup mpl
//
// This specialization of the IfTrue template is selected in case the \a Condition compile time
// constant expression evaluates to \a false. The member type definition is set to the second
// given type \a T2.
*/
template< typename T1    // Type not to be selected
        , typename T2 >  // Type to be selected
struct IfTrue<false,T1,T2>
{
 public:
   //**********************************************************************************************
   typedef T2  Type;  //!< The selected type.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type selection.
// \ingroup mpl
//
// The If class template selects one of the two given types \a T2 and \a T3 depending on \a T1.
// In case \a T1::value evaluates to \a true, the member type definition \a Type is set to \a T2.
// In case \a T1::value evaluates to \a false, \a Type is set to \a T3.
*/
template< typename T1    // Type of the condition
        , typename T2    // Type to be selected if T1::value=true
        , typename T3 >  // Type to be selected if T1::value=false
struct If
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename IfTrue< T1::value, T2, T3 >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
