//=================================================================================================
/*!
//  \file blaze/util/valuetraits/IsPowerOf.h
//  \brief Header file for the IsPowerOf value trait
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

#ifndef _BLAZE_UTIL_VALUETRAITS_ISPOWEROF_H_
#define _BLAZE_UTIL_VALUETRAITS_ISPOWEROF_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for a power relationship of integral values to a given base.
// \ingroup value_traits
//
// This value trait tests whether the given integral value \a N is a power of the base \a B
// according to the equation \f$ B^x = N \f$, where x is any positive integer in the range
// \f$ [0..\infty) \f$. In case the value is a power of \a B, the \a value member enumeration
// is set to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives
// from \a FalseType.

   \code
   blaze::IsPowerOf<2,8>::value   // Evaluates to 1 (2^3 = 8)
   blaze::IsPowerOf<3,27>::value  // Evaluates to 1 (3^3 = 27)
   blaze::IsPowerOf<5,1>::value   // Evaluates to 1 (5^0 = 1)
   blaze::IsPowerOf<1,1>::Type    // Results in TrueType (1^x = 1)
   blaze::IsPowerOf<0,0>          // Is derived from TrueType (0^x = 0)
   blaze::IsPowerOf<2,14>::value  // Evaluates to 0
   blaze::IsPowerOf<1,5>::value   // Evaluates to 0
   blaze::IsPowerOf<0,5>::Type    // Results in FalseType
   blaze::IsPowerOf<2,0>          // Is derived from FalseType
   \endcode
*/
template< size_t B, size_t N >
struct IsPowerOf : public SelectType<N%B,FalseType,typename IsPowerOf<B,N/B>::Type>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = ( N%B )?( 0 ):( IsPowerOf<B,N/B>::value ) };
   typedef typename SelectType<N%B,FalseType,typename IsPowerOf<B,N/B>::Type>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the IsPowerOf value trait.
// \ingroup type_traits
//
// This class ia a partial specialization of the IsPowerOf value trait for any value \a N to
// the base 2. In case \a N is a power of 2, the \a value member enumeration is set to 1, the
// nested type definition \a Type is \a TrueType, and the class derives from \a TrueType. If
// not, \a value is set to 0, \a Type is \a FalseType, and the class derives from \a FalseType.
*/
template< size_t N >
struct IsPowerOf<2,N> : public SelectType<N&(N-1),FalseType,TrueType>::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( N&(N-1) )?( 0 ):( 1 ) };
   typedef typename SelectType<N&(N-1),FalseType,TrueType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the IsPowerOf value trait.
// \ingroup type_traits
//
// This class ia a partial specialization of the IsPowerOf value trait for the value 0 to the
// base 2. Since 0 is no power of 2, this specialization sets the \a value member enumeration
// to 0, the nested type definition \a Type to \a FalseType, and it derives from \a FalseType.
*/
template<>
struct IsPowerOf<2,0> : public FalseType
{
 public:
   //**********************************************************************************************
   enum { value = 0 };
   typedef FalseType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the IsPowerOf value trait.
// \ingroup type_traits
//
// This class ia a partial specialization of the IsPowerOf value trait for the value 1 to any
// given base \a B larger than 1. According to the equation \f$ B^0 = 1 \f$ this specialization
// always sets the \a value member enumeration to 1, the nested type definition \a Type to
// \a TrueType, and it derives from \a TrueType.
*/
template< size_t B >
struct IsPowerOf<B,1> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the IsPowerOf value trait.
// \ingroup type_traits
//
// This class ia a partial specialization of the IsPowerOf value trait for any value \a N larger
// than 1 to the base 1. Since N is no power of 1, this specialization always sets the \a value
// member enumeration to 0, the nested type definition \a Type to \a FalseType, and it derives
// from \a FalseType.
*/
template< size_t N >
struct IsPowerOf<1,N> : public FalseType
{
 public:
   //**********************************************************************************************
   enum { value = 0 };
   typedef FalseType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the IsPowerOf value trait.
// \ingroup type_traits
//
// This class ia a partial specialization of the IsPowerOf value trait for the value 1 to
// the base 1. Since 1 is a power of 1, this specialization always sets the \a value member
// enumeration to 1, the nested type definition \a Type to \a TrueType, and it derives from
// \a TrueType.
*/
template<>
struct IsPowerOf<1,1> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the IsPowerOf value trait.
// \ingroup type_traits
//
// This class ia a partial specialization of the IsPowerOf value trait for the value 0 to the
// base \a B. Since 0 is no power of \a B, this specialization always sets the \a value member
// enumeration to 0, the nested type definition \a Type to \a FalseType, and it derives from
// \a FalseType.
*/
template< size_t B >
struct IsPowerOf<B,0> : public FalseType
{
 public:
   //**********************************************************************************************
   enum { value = 0 };
   typedef FalseType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the IsPowerOf value trait.
// \ingroup type_traits
//
// This class ia a partial specialization of the IsPowerOf value trait for any value \a N to
// the base 0. Since N is no power of 0, this specialization always sets the \a value member
// enumeration to 0, the nested type definition \a Type to \a FalseType, and it derives from
// \a FalseType.
*/
template< size_t N >
struct IsPowerOf<0,N> : public FalseType
{
 public:
   //**********************************************************************************************
   enum { value = 0 };
   typedef FalseType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the IsPowerOf value trait.
// \ingroup type_traits
//
// This class ia a partial specialization of the IsPowerOf value trait for the value 0 to
// the base 0. Since 0 is a power of 0 (\f$ 0^x = 0 \f$), this specialization always sets
// the \a value member enumeration to 1, the nested type definition \a Type to \a TrueType,
// and it derives from \a TrueType.
*/
template<>
struct IsPowerOf<0,0> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
