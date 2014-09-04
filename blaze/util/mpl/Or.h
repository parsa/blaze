//=================================================================================================
/*!
//  \file blaze/util/mpl/Or.h
//  \brief Header file for the Or class template
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_UTIL_MPL_OR_H_
#define _BLAZE_UTIL_MPL_OR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/mpl/Bool.h>
#include <blaze/util/NullType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time logical or evaluation.
// \ingroup mpl
//
// The Or class template performs at compile time a logical or ('&&') evaluation of the up to
// six given compile time conditions:

   \code
   using namespace blaze;

   typedef int  Type;

   Or< IsIntegral<Type>, IsSigned<Type>        >::value  // Evaluates to 1
   Or< IsIntegral<Type>, IsFloatingPoint<Type> >::value  // Evaluates to 1
   Or< IsFloat<Type>   , IsDouble<Type>        >::value  // Evaluates to 0
   \endcode
*/
template< typename T1               // Type of the first operand
        , typename T2               // Type of the second operand
        , typename T3 = NullType    // Type of the third operand
        , typename T4 = NullType    // Type of the fourth operand
        , typename T5 = NullType    // Type of the fifth operand
        , typename T6 = NullType >  // Type of the sixth operand
struct Or
   : public Bool< ( T1::value || T2::value || T3::value || T4::value || T5::value || T6::value ) >
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the Or class template for two operands.
template< typename T1    // Type of the first operand
        , typename T2 >  // Type of the second operand
struct Or<T1,T2,NullType,NullType,NullType,NullType>
   : public Bool< ( T1::value || T2::value ) >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the Or class template for three operands.
template< typename T1    // Type of the first operand
        , typename T2    // Type of the second operand
        , typename T3 >  // Type of the third operand
struct Or<T1,T2,T3,NullType,NullType,NullType>
   : public Bool< ( T1::value || T2::value || T3::value ) >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the Or class template for four operands.
template< typename T1    // Type of the first operand
        , typename T2    // Type of the second operand
        , typename T3    // Type of the third operand
        , typename T4 >  // Type of the fourth operand
struct Or<T1,T2,T3,T4,NullType,NullType>
   : public Bool< ( T1::value || T2::value || T3::value || T4::value ) >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the Or class template for five operands.
template< typename T1    // Type of the first operand
        , typename T2    // Type of the second operand
        , typename T3    // Type of the third operand
        , typename T4    // Type of the fourth operand
        , typename T5 >  // Type of the fifth operand
struct Or<T1,T2,T3,T4,T5,NullType>
   : public Bool< ( T1::value || T2::value || T3::value || T4::value || T5::value ) >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
