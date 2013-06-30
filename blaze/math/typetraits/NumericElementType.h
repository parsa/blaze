//=================================================================================================
/*!
//  \file blaze/math/typetraits/NumericElementType.h
//  \brief Header file for the NumericElementType type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_NUMERICELEMENTTYPE_H_
#define _BLAZE_MATH_TYPETRAITS_NUMERICELEMENTTYPE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the numeric element type of a given data type.
// \ingroup math_type_traits
//
// Via this type trait it is possible to evaluate the numeric (fundamental or complex) element
// type at the heart of a given data type. Examples:

   \code
   blaze::NumericElementType< double >::Type                              // corresponds to double
   blaze::NumericElementType< complex<float> >::Type                      // corresponds to complex<float>
   blaze::NumericElementType< StaticVector<int,3UL> >::Type               // corresponds to int
   blaze::NumericElementType< CompressedVector< DynamicVector<float> > >  // corresponds to float
   \endcode

// Note that per default NumericElementType only supports fundamental/built-in data types,
// complex, and data types with the nested type definition \a ElementType. Support for other
// data types can be added by specializing the NumericElementType class template.
*/
template< typename T >
struct NumericElementType
{
 private:
   //**struct Numeric******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T2 >
   struct Numeric { typedef T2  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct Other********************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T2 >
   struct Other { typedef typename NumericElementType<typename T2::ElementType>::Type  Type; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsNumeric<T>, Numeric<T>, Other<T> >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
