//=================================================================================================
/*!
//  \file blaze/util/typetraits/HasSize.h
//  \brief Header file for the HasSize type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_HASSIZE_H_
#define _BLAZE_UTIL_TYPETRAITS_HASSIZE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time size check.
// \ingroup type_traits
//
// This class offers the possibility to test the size of a type at compile time. If the type
// \a T is exactly \a Size bytes large, the \a value member enumeration is set to 1, the nested
// type definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to 0, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::HasSize<int,4>::value              // Evaluates to 1 (on most architectures)
   blaze::HasSize<float,4>::Type             // Results in TrueType (on most architectures)
   blaze::HasSize<const double,8>            // Is derived from TrueType (on most architectures)
   blaze::HasSize<volatile double,2>::value  // Evaluates to 0
   blaze::HasSize<const char,8>::Type        // Results in FalseType
   blaze::HasSize<unsigned char,4>           // Is derived from FalseType
   \endcode
*/
template< typename T, size_t Size >
class HasSize : public SelectType< sizeof( T ) == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = ( sizeof( T ) == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class ia a partial specialization of the HasSize template for the type \a void. This
// specialization assumes that an object of type \a void has a size of 0. Therefore \a value
// is set to 1, \a Type is \a TrueType, and the class derives from \a TrueType only if the
// \a Size template argument is 0. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< size_t Size >
class HasSize<void,Size> : public SelectType< 0 == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( 0 == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class ia a partial specialization of the HasSize template for constant \a void. This
// specialization assumes that an object of type \a void has a size of 0. Therefore \a value
// is set to 1, \a Type is \a TrueType, and the class derives from \a TrueType only if the
// \a Size template argument is 0. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< size_t Size >
class HasSize<const void,Size> : public SelectType< 0 == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( 0 == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class ia a partial specialization of the HasSize template for volatile \a void. This
// specialization assumes that an object of type \a void has a size of 0. Therefore \a value
// is set to 1, \a Type is \a TrueType, and the class derives from \a TrueType only if the
// \a Size template argument is 0. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< size_t Size >
class HasSize<volatile void,Size> : public SelectType< 0 == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( 0 == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class ia a partial specialization of the HasSize template for constant volatile \a void.
// This specialization assumes that an object of type \a void has a size of 0. Therefore \a value
// is set to 1, \a Type is \a TrueType, and the class derives from \a TrueType only if the \a Size
// template argument is 0. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class
// derives from \a FalseType.
*/
template< size_t Size >
class HasSize<const volatile void,Size> : public SelectType< 0 == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( 0 == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
