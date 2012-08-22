//=================================================================================================
/*!
//  \file blaze/math/traits/TDVecDMatMultTrait.h
//  \brief Header file for the TDVecDMatMultTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_TDVECDMATMULTTRAIT_H_
#define _BLAZE_MATH_TRAITS_TDVECDMATMULTTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Forward.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/SelectType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the expression type of a dense vector/dense matrix multiplication.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the resulting expression type of a dense
// vector/dense matrix multiplication. Given the transpose dense vector type \a VT and the
// row-major dense matrix type \a MT, the nested type \a Type corresponds to the resulting
// expression type. In case either \a VT is not a transpose dense vector type or \a MT is
// not a row-major dense matrix type, the resulting data type \a Type is set to \a INVALID_TYPE.
*/
template< typename VT    // Type of the left-hand side transpose dense vector
        , typename MT >  // Type of the right-hand side row-major dense matrix
struct TDVecDMatMultTrait
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsDenseVector<VT>::value && IsTransposeVector<VT>::value &&
                                IsDenseMatrix<MT>::value && IsRowMajorMatrix<MT>::value
                              , TDVecDMatMultExpr<VT,MT>, INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
