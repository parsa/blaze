//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsResizable.h
//  \brief Header file for the IsResizable type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISRESIZABLE_H_
#define _BLAZE_MATH_TYPETRAITS_ISRESIZABLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/FalseType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for resizable data types.
// \ingroup math_type_traits
//
// This type trait tests whether the given data type is a resizable data type. In case the
// data type can be resized (via the resize() function), the \a value member enumeration
// is set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives
// from \a FalseType. Examples:

   \code
   blaze::IsResizable< DynamicVector<double,false> >::value       // Evaluates to 1
   blaze::IsResizable< const DynamicMatrix<double,false> >::Type  // Results in TrueType
   blaze::IsResizable< volatile CompressedMatrix<int,true> >      // Is derived from TrueType
   blaze::IsResizable< int >::value                               // Evaluates to 0
   blaze::IsResizable< const complex<double> >::Type              // Results in FalseType
   blaze::IsResizable< volatile StaticVector<float,3U,false> >    // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsResizable : public FalseType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = 0 };
   typedef FalseType  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
