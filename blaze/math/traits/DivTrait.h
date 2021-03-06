//=================================================================================================
/*!
//  \file blaze/math/traits/DivTrait.h
//  \brief Header file for the division trait
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_TRAITS_DIVTRAIT_H_
#define _BLAZE_MATH_TRAITS_DIVTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/typetraits/HasDiv.h>
#include <blaze/util/Complex.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/CommonType.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/RemoveCVRef.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, typename, typename = void > struct DivTrait;
template< typename, typename, typename = void > struct DivTraitEval1;
template< typename, typename, typename = void > struct DivTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
auto evalDivTrait( T1&, T2& )
   -> typename DivTraitEval1<T1,T2>::Type;

template< typename T1, typename T2 >
auto evalDivTrait( const T1&, const T2& )
   -> typename DivTrait<T1,T2>::Type;

template< typename T1, typename T2 >
auto evalDivTrait( const volatile T1&, const T2& )
   -> typename DivTrait<T1,T2>::Type;

template< typename T1, typename T2 >
auto evalDivTrait( const T1&, const volatile T2& )
   -> typename DivTrait<T1,T2>::Type;

template< typename T1, typename T2 >
auto evalDivTrait( const volatile T1&, const volatile T2& )
   -> typename DivTrait<T1,T2>::Type;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the DivTrait class.
// \ingroup math_traits
//
// \section divtrait_general General
//
// The DivTrait class template offers the possibility to select the resulting data type of a
// generic division operation between the two given types \a T1 and \a T2. DivTrait defines
// the nested type \a Type, which represents the resulting data type of the division. In case
// the two types \a T1 and \a T2 cannot be divided, a compilation error is created. Note that
// \c const and \c volatile qualifiers and reference modifiers are generally ignored.
//
//
// \n \section divtrait_specializations Creating custom specializations
//
// DivTrait is guaranteed to work for all built-in data types, complex numbers, all vector and
// matrix types of the Blaze library (including views and adaptors) and all data types that
// provide a division operator (i.e. \c operator/). In order to add support for user-defined
// data types that either don't provide a division operator or whose division operator returns
// a proxy object instead of a concrete type (as it is for instance common in expression template
// libraries) it is possible to specialize the DivTrait template. The following example shows the
// according specialization for the division of a dynamic column vector by a double precision
// scalar value:

   \code
   template< typename T1 >
   struct DivTrait< DynamicVector<T1,columnVector>, double >
   {
      using Type = DynamicVector< typename DivTrait<T1,double>::Type, columnVector >;
   };
   \endcode

// \n \section divtrait_examples Examples
//
// The following example demonstrates the use of the DivTrait template, where depending on
// the two given data types the resulting data type is selected:

   \code
   template< typename T1, typename T2 >  // The two generic types
   typename DivTrait<T1,T2>::Type        // The resulting generic return type
   div( const T1& t1, const T2& t2 )     //
   {                                     // The function 'div' returns the
      return t1 / t2;                    // quotient of the two given values
   }                                     //
   \endcode
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct DivTrait
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = decltype( evalDivTrait( std::declval<T1&>(), std::declval<T2&>() ) );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DivTrait class template for two identical builtin types.
// \ingroup math_traits
*/
template< typename T >
struct DivTrait< T, T, EnableIf_t< IsBuiltin_v<T> > >
{
 public:
   //**********************************************************************************************
   using Type = RemoveCVRef_t<T>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DivTrait class template for a complex and a built-in type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct DivTrait< complex<T1>, T2, EnableIf_t< IsBuiltin_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = CommonType_t< complex<T1> , T2 >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DivTrait class template for a built-in and a complex type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct DivTrait< T1, complex<T2>, EnableIf_t< IsBuiltin_v<T1> > >
{
 public:
   //**********************************************************************************************
   using Type = CommonType_t< T1, complex<T2> >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DivTrait class template for two complex types.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct DivTrait< complex<T1>, complex<T2> >
{
 public:
   //**********************************************************************************************
   using Type = CommonType_t< complex<T1>, complex<T2> >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the DivTrait class template.
// \ingroup math_traits
//
// The DivTrait_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the DivTrait class template. For instance, given the types \a T1 and \a T2 the following
// two type definitions are identical:

   \code
   using Type1 = typename blaze::DivTrait<T1,T2>::Type;
   using Type2 = blaze::DivTrait_t<T1,T2>;
   \endcode
*/
template< typename T1, typename T2 >
using DivTrait_t = typename DivTrait<T1,T2>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the DivTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct DivTraitEval1
{
 public:
   //**********************************************************************************************
   using Type = typename DivTraitEval2<T1,T2>::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the DivTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct DivTraitEval2
{
 private:
   //**********************************************************************************************
   struct DivType { using Type = decltype( std::declval<T1>() / std::declval<T2>() ); };
   struct Failure { using Type = INVALID_TYPE; };
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   using Type = typename If_t< HasDiv_v<T1,T2>, DivType, Failure >::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
