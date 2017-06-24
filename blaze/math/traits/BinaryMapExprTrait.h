//=================================================================================================
/*!
//  \file blaze/math/traits/BinaryMapExprTrait.h
//  \brief Header file for the BinaryMapExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_BINARYMAPEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_BINARYMAPEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/DMatDMatMapExprTrait.h>
#include <blaze/math/traits/DMatTDMatMapExprTrait.h>
#include <blaze/math/traits/DVecDVecMapExprTrait.h>
#include <blaze/math/traits/TDMatDMatMapExprTrait.h>
#include <blaze/math/traits/TDMatTDMatMapExprTrait.h>
#include <blaze/math/traits/TDVecTDVecMapExprTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/util/InvalidType.h>
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
/*!\brief Evaluation of the return type of a binary map expression.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a binary map expression.
// Given the two types \a T1 and \a T2, which must either be vector or matrix types, and the
// custom operation type \a OP, the nested type \a Type corresponds to the resulting return
// type. In case the types of \a T1 or \a T2 don't fit or if no binary map operation exists
// for the types, the resulting data type \a Type is set to \a INVALID_TYPE.
*/
template< typename T1    // Type of the left-hand side binary map operand
        , typename T2    // Type of the right-hand side binary map operand
        , typename OP >  // Type of the custom operation
struct BinaryMapExprTrait
{
 private:
   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { using Type = INVALID_TYPE; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Tmp = If_< IsDenseMatrix<T1>
                  , If_< IsDenseMatrix<T2>
                       , If_< IsRowMajorMatrix<T1>
                            , If_< IsRowMajorMatrix<T2>
                                 , DMatDMatMapExprTrait<T1,T2,OP>
                                 , DMatTDMatMapExprTrait<T1,T2,OP> >
                            , If_< IsRowMajorMatrix<T2>
                                 , TDMatDMatMapExprTrait<T1,T2,OP>
                                 , TDMatTDMatMapExprTrait<T1,T2,OP> > >
                       , Failure >
                  , If_< IsDenseVector<T1>
                       , If_< IsDenseVector<T2>
                            , If_< IsRowVector<T1>
                                 , If_< IsRowVector<T2>
                                      , TDVecTDVecMapExprTrait<T1,T2,OP>
                                      , Failure >
                                 , If_< IsRowVector<T2>
                                      , Failure
                                      , DVecDVecMapExprTrait<T1,T2,OP> > >
                            , Failure >
                       , Failure >
                  >;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<T1>, IsVolatile<T1>, IsReference<T1>
                                , IsConst<T2>, IsVolatile<T2>, IsReference<T2> >
                            , BinaryMapExprTrait< Decay_<T1>, Decay_<T2>, OP >
                            , Tmp >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the BinaryMapExprTrait class template.
// \ingroup math_traits
//
// The BinaryMapExprTrait_ alias declaration provides a convenient shortcut to access the nested
// \a Type of the BinaryMapExprTrait class template. For instance, given the data types \a T1 and
// and \a T2 and the custom operation type \a OP the following two type definitions are identical:

   \code
   using Type1 = typename BinaryMapExprTrait<T1,T2,OP>::Type;
   using Type2 = BinaryMapExprTrait_<T1,T2,OP>;
   \endcode
*/
template< typename T1    // Type of the left-hand side map operand
        , typename T2    // Type of the right-hand side map operand
        , typename OP >  // Type of the custom operation
using BinaryMapExprTrait_ = typename BinaryMapExprTrait<T1,T2,OP>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
