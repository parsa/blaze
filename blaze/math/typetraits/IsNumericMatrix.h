//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsNumericMatrix.h
//  \brief Header file for the IsNumericMatrix type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISNUMERICMATRIX_H_
#define _BLAZE_MATH_TYPETRAITS_ISNUMERICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for numeric matrix types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a numeric matrix type,
// i.e. a matrix with numeric element type. In case the type is a numeric matrix type, the
// \a value member constant is set to \a true, the nested type definition \a Type is \a TrueType,
// and the class derives from \a TrueType. Otherwise \a yes is set to \a false, \a Type is
// \a FalseType, and the class derives from \a FalseType.

   \code
   using Type1 = DynamicMatrix<int>;
   using Type2 = CompressedMatrix< complex<double> >;
   using Type3 = LowerMatrix< DynamicMatrix<float> >;

   using Type4 = double;
   using Type5 = DynamicVector<int>;
   using Type6 = DynamicMatrix< DynamicVector<int> >;

   blaze::IsNumericMatrix< Type1 >::value  // Evaluates to 1
   blaze::IsNumericMatrix< Type2 >::Type   // Results in TrueType
   blaze::IsNumericMatrix< Type3 >         // Is derived from TrueType
   blaze::IsNumericMatrix< Type4 >::value  // Evaluates to 0
   blaze::IsNumericMatrix< Type5 >::Type   // Results in FalseType
   blaze::IsNumericMatrix< Type6 >         // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsNumericMatrix
   : public BoolConstant< IsMatrix_v<T> && IsNumeric_v< UnderlyingElement_t<T> > >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsNumericMatrix type trait.
// \ingroup math_type_traits
//
// The IsNumericMatrix_v variable template provides a convenient shortcut to access the nested
// \a value of the IsNumericMatrix class template. For instance, given the type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = blaze::IsNumericMatrix<T>::value;
   constexpr bool value2 = blaze::IsNumericMatrix_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsNumericMatrix_v = IsNumericMatrix<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
