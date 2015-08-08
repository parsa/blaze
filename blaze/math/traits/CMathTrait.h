//=================================================================================================
/*!
//  \file blaze/math/traits/CMathTrait.h
//  \brief Header file for the cmath trait
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

#ifndef _BLAZE_MATH_TRAITS_CMATHTRAIT_H_
#define _BLAZE_MATH_TRAITS_CMATHTRAIT_H_


namespace blaze {

//=================================================================================================
//
//  CMATH TRAIT
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the CMathTrait class.
// \ingroup math_traits
//
// The CMathTrait template evaluates the return type of the mathematical functions defined
// in the C++ header \<cmath\> depending on the type of the template argument. In case of an
// integral data type or double precision argument, the return value of the functions is
// double, whereas the return type is float for single precision arguments and long double
// for long double precision arguments.
//
// <table border="0" cellspacing="0" cellpadding="1">
//    <tr>
//       <td width="250px"> <b> Template argument T </b> </td>
//       <td width="100px"> <b> Type </b> </td>
//    </tr>
//    <tr>
//       <td>float</td>
//       <td>float</td>
//    </tr>
//    <tr>
//       <td>integral data types and double</td>
//       <td>double</td>
//    </tr>
//    <tr>
//       <td>long double</td>
//       <td>long double</td>
//    </tr>
// </table>
*/
template< typename T >
struct CMathTrait
{
   typedef double Type;  //!< Return type of the \<cmath\> functions for integral and double arguments.
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief CMathTrait<float> specialization.
// \ingroup math_traits
*/
template<>
struct CMathTrait<float>
{
   typedef float Type;  //!< Return type of the \<cmath\> functions for float arguments.
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief CMathTrait<long double> specialization.
// \ingroup math_traits
*/
template<>
struct CMathTrait<long double>
{
   typedef long double Type;  //!< Return type of the \<cmath\> functions for long double arguments.
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
