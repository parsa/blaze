//=================================================================================================
/*!
//  \file blaze/math/traits/BinaryMapTrait.h
//  \brief Header file for the binary map trait
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_TRAITS_BINARYMAPTRAIT_H_
#define _BLAZE_MATH_TRAITS_BINARYMAPTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/util/typetraits/Decay.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, typename, typename, typename = void > struct BinaryMapTrait;
template< typename, typename, typename, typename = void > struct BinaryMapTraitEval1;
template< typename, typename, typename, typename = void > struct BinaryMapTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
auto evalBinaryMapTrait( T1&, T2&, OP )
   -> typename BinaryMapTraitEval1<T1,T2,OP>::Type;

template< typename T1, typename T2, typename OP >
auto evalBinaryMapTrait( const T1&, const T2&, OP )
   -> typename BinaryMapTrait<T1,T2,OP>::Type;

template< typename T1, typename T2, typename OP >
auto evalBinaryMapTrait( const T1&, volatile const T2&, OP )
   -> typename BinaryMapTrait<T1,T2,OP>::Type;

template< typename T1, typename T2, typename OP >
auto evalBinaryMapTrait( volatile const T1&, const T2&, OP )
   -> typename BinaryMapTrait<T1,T2,OP>::Type;

template< typename T1, typename T2, typename OP >
auto evalBinaryMapTrait( volatile const T1&, volatile const T2&, OP )
   -> typename BinaryMapTrait<T1,T2,OP>::Type;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the BinaryMapTrait class.
// \ingroup math_traits
//
// \section binarymaptrait_general General
//
// The BinaryMapTrait class template offers the possibility to select the resulting data type of
// a generic, binary map operation between the two given types \a T1 and \a T2. BinaryMapTrait
// defines the nested type \a Type, which represents the resulting data type of the map operation.
// In case no result type can be determined for the two types \a T1 and \a T2, a compilation error
// is created. Note that \c const and \c volatile qualifiers and reference modifiers are generally
// ignored.
//
//
// \n \section binarymaptrait_specializations Creating custom specializations
//
// BinaryMapTrait is guaranteed to work for all built-in data types, complex numbers, all vector
// and matrix types of the Blaze library (including views and adaptors) and all data types that
// work in combination with the provided custom operation \a OP. In order to add support for
// user-defined data types or in order to adapt to special cases it is possible to specialize
// the BinaryMapTrait template. The following example shows the according specialization for
// binary map operations between two dynamic column vectors:

   \code
   template< typename T1, typename T2, typename OP >
   struct BinaryMapTrait< DynamicVector<T1,columnVector>, DynamicVector<T2,columnVector>, OP >
   {
      using Type = DynamicVector< typename BinaryMapTrait<T1,T2,OP>::Type, columnVector >;
   };
   \endcode
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct BinaryMapTrait
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = decltype( evalBinaryMapTrait( std::declval<T1&>(), std::declval<T2&>(), std::declval<OP&>() ) );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the BinaryMapTrait class template.
// \ingroup math_traits
//
// The BinaryMapTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the BinaryMapTrait class template. For instance, given the types \a T1 and \a T2
// and the custom operation type \a OP the following two type definitions are identical:

   \code
   using Type1 = typename blaze::BinaryMapTrait<T1,T2,OP>::Type;
   using Type2 = blaze::BinaryMapTrait_t<T1,T2,OP>;
   \endcode
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2    // Type of the right-hand side operand
        , typename OP >  // Type of the custom operation
using BinaryMapTrait_t = typename BinaryMapTrait<T1,T2,OP>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the BinaryMapTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct BinaryMapTraitEval1
{
 public:
   //**********************************************************************************************
   using Type = typename BinaryMapTraitEval2<T1,T2,OP>::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the BinaryMapTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct BinaryMapTraitEval2
{
 public:
   //**********************************************************************************************
   using Type = Decay_t< decltype( std::declval<OP>()( std::declval<T1>(), std::declval<T2>() ) ) >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
