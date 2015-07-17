//=================================================================================================
/*!
//  \file blaze/math/traits/ConjExprTrait.h
//  \brief Header file for the ConjExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_CONJEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_CONJEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/DMatConjExprTrait.h>
#include <blaze/math/traits/DVecConjExprTrait.h>
#include <blaze/math/traits/SMatConjExprTrait.h>
#include <blaze/math/traits/SVecConjExprTrait.h>
#include <blaze/math/traits/TDMatConjExprTrait.h>
#include <blaze/math/traits/TDVecConjExprTrait.h>
#include <blaze/math/traits/TSMatConjExprTrait.h>
#include <blaze/math/traits/TSVecConjExprTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>
#include <blaze/util/typetraits/RemoveCV.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the return type of a complex conjugate expression.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a complex conjugate
// expression. Given the type \a T, which must either be a scalar, vector, or matrix type,
// the nested type \a Type corresponds to the resulting return type. In case the type of
// \a T doesn't fit or if no complex conjugate operation exists for the type, the resulting
// data type \a Type is set to \a INVALID_TYPE.
*/
template< typename T >  // Type of the complex conjugate operand
struct ConjExprTrait
{
 private:
   //**struct ScalarConj***************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename ST >
   struct ScalarConj { typedef T  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsMatrix<T>
                      , typename If< IsDenseMatrix<T>
                                   , typename If< IsRowMajorMatrix<T>
                                                , DMatConjExprTrait<T>
                                                , TDMatConjExprTrait<T>
                                                >::Type
                                   , typename If< IsRowMajorMatrix<T>
                                                , SMatConjExprTrait<T>
                                                , TSMatConjExprTrait<T>
                                                >::Type
                                   >::Type
                      , typename If< IsVector<T>
                                   , typename If< IsDenseVector<T>
                                                , typename If< IsRowVector<T>
                                                             , TDVecConjExprTrait<T>
                                                             , DVecConjExprTrait<T>
                                                             >::Type
                                                , typename If< IsRowVector<T>
                                                             , TSVecConjExprTrait<T>
                                                             , SVecConjExprTrait<T>
                                                             >::Type
                                                >::Type
                                   , typename If< IsNumeric<T>
                                                , ScalarConj<T>
                                                , Failure
                                                >::Type
                                   >::Type
                      >::Type  Tmp;

   typedef typename RemoveReference< typename RemoveCV<T>::Type >::Type  Type1;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< Or< IsConst<T>, IsVolatile<T>, IsReference<T> >
                      , ConjExprTrait<Type1>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
