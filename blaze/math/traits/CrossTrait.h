//=================================================================================================
/*!
//  \file blaze/math/traits/CrossTrait.h
//  \brief Header file for the cross product trait
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

#ifndef _BLAZE_MATH_TRAITS_CROSSTRAIT_H_
#define _BLAZE_MATH_TRAITS_CROSSTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/InvalidType.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsVolatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the CrossTrait class.
// \ingroup math
//
// \section crosstrait_general General
//
// The CrossTrait class template offers the possibility to select the resulting data type of
// a generic cross product operation between the two given types \a T1 and \a T2. CrossTrait
// defines the nested type \a Type, which represents the resulting data type of the cross
// product. In case \a T1 and \a T2 cannot be combined in a cross product, the resulting data
// type \a Type is set to \a INVALID_TYPE. Note that \a const and \a volatile qualifiers and
// reference modifiers are generally ignored.
//
// Since the cross product is only defined for 3-dimensional vectors, the CrossTrait template
// only supports the following vector types:
//
// <ul>
//    <li>blaze::StaticVector</li>
//    <li>blaze::DynamicVector</li>
//    <li>blaze::CompressedVector</li>
// </ul>
//
//
// \n \section specializations Creating custom specializations
//
// It is possible to specialize the CrossTrait template for additional user-defined data types.
// The following example shows the according specialization for the cross product between two
// static column vectors:

   \code
   template< typename T1, typename T2 >
   struct CrossTrait< StaticVector<T1,3UL,false>, StaticVector<T2,3UL,false> >
   {
      typedef StaticVector< typename SubTrait< typename MultTrait<T1,T2>::Type
                                             , typename MultTrait<T1,T2>::Type >::Type, 3UL, false >  Type;
   };
   \endcode

// \n \section crosstrait_examples Examples
//
// The following example demonstrates the use of the CrossTrait template, where depending on
// the two given data types the resulting data type is selected:

   \code
   template< typename T1, typename T2 >    // The two generic types
   typename CrossTrait<T1,T2>::Type        // The resulting generic return type
   cross( T1 t1, T2 t2 )                   //
   {                                       // The function 'cross' returns the cross
      return t1 % t2;                      // product of the two given values
   }                                       //
   \endcode
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
struct CrossTrait
{
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef INVALID_TYPE  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
