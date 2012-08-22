//=================================================================================================
/*!
//  \file blaze/math/CMathTrait.h
//  \brief Header file for the cmath trait
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

#ifndef _BLAZE_MATH_CMATHTRAIT_H_
#define _BLAZE_MATH_CMATHTRAIT_H_


namespace blaze {

//=================================================================================================
//
//  CMATH TRAIT
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the CMathTrait class.
// \ingroup math
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
// \ingroup math
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
// \ingroup math
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
