//=================================================================================================
/*!
//  \file blaze/math/traits/ImagExprTrait.h
//  \brief Header file for the ImagExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_IMAGEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_IMAGEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/DMatImagExprTrait.h>
#include <blaze/math/traits/DVecImagExprTrait.h>
#include <blaze/math/traits/SMatImagExprTrait.h>
#include <blaze/math/traits/SVecImagExprTrait.h>
#include <blaze/math/traits/TDMatImagExprTrait.h>
#include <blaze/math/traits/TDVecImagExprTrait.h>
#include <blaze/math/traits/TSMatImagExprTrait.h>
#include <blaze/math/traits/TSVecImagExprTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/Decay.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the return type of a imaginary part expression.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a imaginary part expression.
// Given the type \a T, which must either be a scalar, vector, or matrix type, the nested type
// \a Type corresponds to the resulting return type. In case the type of \a T doesn't fit or
// if no imaginary part operation exists for the type, the resulting data type \a Type is set
// to \a INVALID_TYPE.
*/
template< typename T >  // Type of the imaginary part operand
struct ImagExprTrait
{
 private:
   //**struct Builtin******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename ST >
   struct Builtin { using Type = T; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct Complex******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename CT >
   struct Complex { using Type = typename CT::value_type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { using Type = INVALID_TYPE; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Tmp = If_< IsMatrix<T>
                  , If_< IsDenseMatrix<T>
                       , If_< IsRowMajorMatrix<T>
                            , DMatImagExprTrait<T>
                            , TDMatImagExprTrait<T> >
                       , If_< IsRowMajorMatrix<T>
                            , SMatImagExprTrait<T>
                            , TSMatImagExprTrait<T> > >
                  , If_< IsVector<T>
                       , If_< IsDenseVector<T>
                            , If_< IsRowVector<T>
                                 , TDVecImagExprTrait<T>
                                 , DVecImagExprTrait<T> >
                            , If_< IsRowVector<T>
                                 , TSVecImagExprTrait<T>
                                 , SVecImagExprTrait<T> > >
                       , If_< IsNumeric<T>
                            , If_< IsComplex<T>
                                 , Complex<T>
                                 , Builtin<T> >
                            , Failure > > >;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<T>, IsVolatile<T>, IsReference<T> >
                            , ImagExprTrait< Decay_<T> >
                            , Tmp >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the ImagExprTrait class template.
// \ingroup math_traits
//
// The ImagExprTrait_ alias declaration provides a convenient shortcut to access the nested
// \a Type of the ImagExprTrait class template. For instance, given the type \a T the following
// two type definitions are identical:

   \code
   using Type1 = typename ImagExprTrait<T>::Type;
   using Type2 = ImagExprTrait_<T>;
   \endcode
*/
template< typename T >  // Type of the imaginary part operand
using ImagExprTrait_ = typename ImagExprTrait<T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
