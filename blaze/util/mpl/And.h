//=================================================================================================
/*!
//  \file blaze/util/mpl/And.h
//  \brief Header file for the And class template
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

#ifndef _BLAZE_UTIL_MPL_AND_H_
#define _BLAZE_UTIL_MPL_AND_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/NullType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time logical and evaluation.
// \ingroup mpl
//
// The And class template performs at compile time a logical and ('&&') evaluation of the up to
// five given compile time conditions:

   \code
   using namespace blaze;

   typedef int  Type;

   And< IsIntegral<Type>, IsSigned<Type>        >::value  // Evaluates to 1
   And< IsIntegral<Type>, IsFloatingPoint<Type> >::value  // Evaluates to 0
   And< IsFloat<Type>   , IsDouble<Type>        >::value  // Evaluates to 0
   \endcode
*/
template< typename T1               // Type of the first operand
        , typename T2               // Type of the second operand
        , typename T3 = NullType    // Type of the third operand
        , typename T4 = NullType    // Type of the fourth operand
        , typename T5 = NullType >  // Type of the fifth operand
struct And
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = T1::value && T2::value && T3::value && T4::value && T5::value };
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the And class template for two operands.
template< typename T1    // Type of the first operand
        , typename T2 >  // Type of the second operand
struct And<T1,T2,NullType,NullType,NullType>
{
 public:
   //**********************************************************************************************
   enum { value = T1::value && T2::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the And class template for three operands.
template< typename T1    // Type of the first operand
        , typename T2    // Type of the second operand
        , typename T3 >  // Type of the third operand
struct And<T1,T2,T3,NullType,NullType>
{
 public:
   //**********************************************************************************************
   enum { value = T1::value && T2::value && T3::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the And class template for four operands.
template< typename T1    // Type of the first operand
        , typename T2    // Type of the second operand
        , typename T3    // Type of the third operand
        , typename T4 >  // Type of the fourth operand
struct And<T1,T2,T3,T4,NullType>
{
 public:
   //**********************************************************************************************
   enum { value = T1::value && T2::value && T3::value && T4::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
