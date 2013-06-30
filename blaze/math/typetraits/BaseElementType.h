//=================================================================================================
/*!
//  \file blaze/math/typetraits/BaseElementType.h
//  \brief Header file for the BaseElementType type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_BASEELEMENTTYPE_H_
#define _BLAZE_MATH_TYPETRAITS_BASEELEMENTTYPE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the base element type of a given data type.
// \ingroup math_type_traits
//
// Via this type trait it is possible to evaluate the base (fundamental) element type at the
// heart of a given data type. Examples:

   \code
   blaze::BaseElementType< double >::Type                              // corresponds to double
   blaze::BaseElementType< complex<float> >::Type                      // corresponds to float
   blaze::BaseElementType< StaticVector<int,3UL> >::Type               // corresponds to int
   blaze::BaseElementType< CompressedVector< DynamicVector<float> > >  // corresponds to float
   \endcode

// Note that per default BaseElementType only supports fundamental/built-in data types, complex,
// and data types with the nested type definition \a ElementType. Support for other data types
// can be added by specializing the BaseElementType class template.
*/
template< typename T >
struct BaseElementType
{
 private:
   //**struct Builtin******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T2 >
   struct Builtin { typedef T2  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct Complex******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T2 >
   struct Complex { typedef typename BaseElementType<typename T2::value_type>::Type  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct Other********************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T2 >
   struct Other { typedef typename BaseElementType<typename T2::ElementType>::Type  Type; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsBuiltin<T>
                     , Builtin<T>
                     , typename If< IsComplex<T>
                                  , Complex<T>
                                  , Other<T>
                                  >::Type
                     >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
