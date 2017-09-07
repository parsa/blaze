//=================================================================================================
/*!
//  \file blaze/math/traits/UnaryMapTrait.h
//  \brief Header file for the unary map trait
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_TRAITS_UNARYMAPTRAIT_H_
#define _BLAZE_MATH_TRAITS_UNARYMAPTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
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
/*!\brief Base template for the UnaryMapTrait class.
// \ingroup math_traits
//
// \section unarymaptrait_general General
//
// The UnaryMapTrait class template offers the possibility to select the resulting data type of
// a generic, unary map operation on the given type \a T. UnaryMapTrait defines the nested type
// \a Type, which represents the resulting data type of the map operation. In case no result type
// can be determined for the type \a T, a compilation error is created. Note that \c const and
// \c volatile qualifiers and reference modifiers are generally ignored.
//
//
// \n \section unarymaptrait_specializations Creating custom specializations
//
// UnaryMapTrait is guaranteed to work for all built-in data types, complex numbers, all vector
// and matrix types of the Blaze library (including views and adaptors) and all data types that
// work in combination with the provided custom operation \a OP. In order to add support for
// user-defined data types or in order to adapt to special cases it is possible to specialize
// the UnaryMapTrait template. The following example shows the according specialization for
// unary map operations with a dynamic column vector:

   \code
   template< typename T, typename OP >
   struct UnaryMapTrait< DynamicVector<T,columnVector>, OP >
   {
      using Type = DynamicVector< typename UnaryMapTrait<T,OP>::Type, columnVector >;
   };
   \endcode
*/
template< typename T     // Type of the operand
        , typename OP >  // Type of the custom operation
struct UnaryMapTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct MappedType
   {
      using Type = Decay_< decltype( std::declval<OP>()( std::declval<T>() ) ) >;
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<T>, IsVolatile<T>, IsReference<T>, IsAdaptor<T> >
                            , UnaryMapTrait< RemoveAdaptor_< Decay_<T> >, OP >
                            , MappedType >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the UnaryMapTrait class template.
// \ingroup math_traits
//
// The UnaryMapTrait_ alias declaration provides a convenient shortcut to access the nested \a Type
// of the UnaryMapTrait class template. For instance, given the type \a T and the custom operation
// type \a OP the following two type definitions are identical:

   \code
   using Type1 = typename UnaryMapTrait<T,OP>::Type;
   using Type2 = UnaryMapTrait_<T,OP>;
   \endcode
*/
template< typename T     // Type of the operand
        , typename OP >  // Type of the custom operation
using UnaryMapTrait_ = typename UnaryMapTrait<T,OP>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
