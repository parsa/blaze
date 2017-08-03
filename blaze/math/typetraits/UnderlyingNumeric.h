//=================================================================================================
/*!
//  \file blaze/math/typetraits/UnderlyingNumeric.h
//  \brief Header file for the UnderlyingNumeric type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_UNDERLYINGNUMERIC_H_
#define _BLAZE_MATH_TYPETRAITS_UNDERLYINGNUMERIC_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the underlying numeric element type of a given data type.
// \ingroup math_type_traits
//
// Via this type trait it is possible to evaluate the underlying numeric (fundamental or complex)
// element type at the heart of a given data type. Examples:

   \code
   using Type1 = double;                                    // Built-in data type
   using Type2 = complex<float>;                            // Complex data type
   using Type3 = StaticVector<int,3UL>;                     // Vector with built-in element type
   using Type4 = CompressedVector< DynamicVector<float> >;  // Vector with vector element type

   blaze::UnderlyingNumeric< Type1 >::Type  // corresponds to double
   blaze::UnderlyingNumeric< Type2 >::Type  // corresponds to complex<float>
   blaze::UnderlyingNumeric< Type3 >::Type  // corresponds to int
   blaze::UnderlyingNumeric< Type4 >::Type  // corresponds to float
   \endcode

// Note that per default UnderlyingNumeric only supports fundamental/built-in data types, complex,
// and data types with the nested type definition \a ElementType. Support for other data types can
// be added by specializing the UnderlyingNumeric class template.
*/
template< typename T >
struct UnderlyingNumeric
{
 private:
   //**struct BuiltinOrComplex*********************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T2 >
   struct BuiltinOrComplex { using Type = T2; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct Other********************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T2 >
   struct Other { using Type = typename UnderlyingNumeric<typename T2::ElementType>::Type; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsBuiltin<T>, IsComplex<T> >
                            , BuiltinOrComplex<T>
                            , Other<T>
                            >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the UnderlyingNumeric type trait.
// \ingroup type_traits
//
// The UnderlyingNumeric_ alias declaration provides a convenient shortcut to access the
// nested \a Type of the UnderlyingNumeric class template. For instance, given the type \a T
// the following two type definitions are identical:

   \code
   using Type1 = typename UnderlyingNumeric<T>::Type;
   using Type2 = UnderlyingNumeric_<T>;
   \endcode
*/
template< typename T >
using UnderlyingNumeric_ = typename UnderlyingNumeric<T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
