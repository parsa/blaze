//=================================================================================================
/*!
//  \file blaze/math/traits/CrossExprTrait.h
//  \brief Header file for the CrossExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_CROSSEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_CROSSEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/DVecDVecCrossExprTrait.h>
#include <blaze/math/traits/DVecSVecCrossExprTrait.h>
#include <blaze/math/traits/SVecDVecCrossExprTrait.h>
#include <blaze/math/traits/SVecSVecCrossExprTrait.h>
#include <blaze/math/typetraits/IsColumnVector.h>
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
/*!\brief Evaluation of the return type of a cross product expression.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a cross product expression.
// Given the two types \a T1 and \a T2, which must be dense or sparse column vectors, the nested
// type \a Type corresponds to the resulting return type. In case the types of \a T1 or \a T2
// don't fit or if the two types cannot be used in a cross product, the resulting data type
// \a Type is set to \a INVALID_TYPE.
*/
template< typename T1    // Type of the left-hand side cross product operand
        , typename T2 >  // Type of the right-hand side cross product operand
struct CrossExprTrait
{
 private:
   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsVector<T1>
                      , typename If< IsVector<T2>
                                   , typename If< IsColumnVector<T1>
                                                , typename If< IsColumnVector<T2>
                                                             , typename If< IsDenseVector<T1>
                                                                          , typename If< IsDenseVector<T2>
                                                                                       , DVecDVecCrossExprTrait<T1,T2>
                                                                                       , DVecSVecCrossExprTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsDenseVector<T2>
                                                                                       , SVecDVecCrossExprTrait<T1,T2>
                                                                                       , SVecSVecCrossExprTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , Failure
                                                             >::Type
                                                , Failure
                                                >::Type
                                   , Failure
                                   >::Type
                      , Failure
                      >::Type  Tmp;

   typedef typename RemoveReference< typename RemoveCV<T1>::Type >::Type  Type1;
   typedef typename RemoveReference< typename RemoveCV<T2>::Type >::Type  Type2;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< Or< IsConst<T1>, IsVolatile<T1>, IsReference<T1>
                          , IsConst<T2>, IsVolatile<T2>, IsReference<T2> >
                      , CrossExprTrait<Type1,Type2>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
