//=================================================================================================
/*!
//  \file blaze/math/traits/BinaryMapTrait.h
//  \brief Header file for the binary map trait
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

#ifndef _BLAZE_MATH_TRAITS_BINARYMAPTRAIT_H_
#define _BLAZE_MATH_TRAITS_BINARYMAPTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/Decay.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

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
// Per default, BinaryMapTrait supports all built-in data types. Additionally, the Blaze library
// provides appropriate specializations for the following user-defined arithmetic types:
//
// <ul>
//    <li>std::complex</li>
//    <li>blaze::StaticVector</li>
//    <li>blaze::HybridVector</li>
//    <li>blaze::DynamicVector</li>
//    <li>blaze::CustomVector</li>
//    <li>blaze::StaticMatrix</li>
//    <li>blaze::HybridMatrix</li>
//    <li>blaze::DynamicMatrix</li>
//    <li>blaze::CustomMatrix</li>
//    <li>blaze::SymmetricMatrix</li>
//    <li>blaze::HermitianMatrix</li>
//    <li>blaze::LowerMatrix</li>
//    <li>blaze::UniLowerMatrix</li>
//    <li>blaze::StrictlyLowerMatrix</li>
//    <li>blaze::UpperMatrix</li>
//    <li>blaze::UniUpperMatrix</li>
//    <li>blaze::StrictlyUpperMatrix</li>
//    <li>blaze::DiagonalMatrix</li>
// </ul>
//
//
// \n \section binarymaptrait_specializations Creating custom specializations
//
// BinaryMapTrait is guaranteed to work for all data types that work in combination with the
// provided custom operation \a OP. In order to add support for user-defined data types or in
// order to adapt to special cases it is possible to specialize the BinaryMapTrait template.
// The following example shows the according specialization for binary map operations between
// two dynamic column vectors:

   \code
   template< typename T1, typename T2, typename OP >
   struct BinaryMapTrait< DynamicVector<T1,columnVector>, DynamicVector<T2,columnVector>, OP >
   {
      using Type = DynamicVector< typename BinaryMapTrait<T1,T2,OP>::Type, columnVector >;
   };
   \endcode
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2    // Type of the right-hand side operand
        , typename OP >  // Type of the custom operation
struct BinaryMapTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct MappedType
   {
      using Type = Decay_< decltype( std::declval<OP>()( std::declval<T1>(), std::declval<T2>() ) ) >;
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<T1>, IsVolatile<T1>, IsReference<T1>, IsAdaptor<T1>
                                , IsConst<T2>, IsVolatile<T2>, IsReference<T2>, IsAdaptor<T2> >
                            , BinaryMapTrait< RemoveAdaptor_< Decay_<T1> >
                                            , RemoveAdaptor_< Decay_<T2> >, OP >
                            , MappedType >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the BinaryMapTrait class template.
// \ingroup math_traits
//
// The BinaryMapTrait_ alias declaration provides a convenient shortcut to access the nested
// \a Type of the BinaryMapTrait class template. For instance, given the types \a T1 and \a T2
// and the custom operation type \a OP the following two type definitions are identical:

   \code
   using Type1 = typename BinaryMapTrait<T1,T2,OP>::Type;
   using Type2 = BinaryMapTrait_<T1,T2,OP>;
   \endcode
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2    // Type of the right-hand side operand
        , typename OP >  // Type of the custom operation
using BinaryMapTrait_ = typename BinaryMapTrait<T1,T2,OP>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
