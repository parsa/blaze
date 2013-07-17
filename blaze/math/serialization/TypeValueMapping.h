//=================================================================================================
/*!
//  \file blaze/math/serialization/TypeValueMapping.h
//  \brief Header file for the TypeValueMapping class template
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

#ifndef _BLAZE_MATH_SERIALIZATION_TypeValueMapping_H_
#define _BLAZE_MATH_SERIALIZATION_TypeValueMapping_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/IsSigned.h>
#include <blaze/util/typetraits/IsUnsigned.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the TypeValueMapping class template.
// \ingroup math_serialization
*/
template< bool IsSignedIntegral, bool IsUnsignedIntegral, bool IsFloatingPoint, bool IsComplex >
struct TypeValueMappingHelper;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TypeValueMappingHelper for compound data types.
// \ingroup math_serialization
*/
template<>
struct TypeValueMappingHelper<false,false,false,false>
{
 public:
   //**********************************************************************************************
   enum { value = 0 };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TypeValueMappingHelper for signed integral data types.
// \ingroup math_serialization
*/
template<>
struct TypeValueMappingHelper<true,false,false,false>
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TypeValueMappingHelper for unsigned integral data types.
// \ingroup math_serialization
*/
template<>
struct TypeValueMappingHelper<false,true,false,false>
{
 public:
   //**********************************************************************************************
   enum { value = 2 };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TypeValueMappingHelper for floating-point data types.
// \ingroup math_serialization
*/
template<>
struct TypeValueMappingHelper<false,false,true,false>
{
 public:
   //**********************************************************************************************
   enum { value = 3 };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TypeValueMappingHelper for complex data types.
// \ingroup math_serialization
*/
template<>
struct TypeValueMappingHelper<false,false,false,true>
{
 public:
   //**********************************************************************************************
   enum { value = 4 };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion from a data type to a serial representation.
// \ingroup math_serialization
//
// This class template converts the given data type into an integral representation suited for
// serialization. Depending on the given data type, the \a value member enumeration is set to
// the according serial representation.
*/
template< typename T >
struct TypeValueMapping
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = TypeValueMappingHelper< IsIntegral<T>::value && IsSigned<T>::value
                                        , IsIntegral<T>::value && IsUnsigned<T>::value
                                        , IsFloatingPoint<T>::value
                                        , IsComplex<T>::value
                                        >::value };
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
