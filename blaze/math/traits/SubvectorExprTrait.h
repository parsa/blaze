//=================================================================================================
/*!
//  \file blaze/math/traits/SubvectorExprTrait.h
//  \brief Header file for the SubvectorExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_SUBVECTOREXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_SUBVECTOREXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsTransExpr.h>
#include <blaze/math/views/Forward.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/IsConst.h>
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
/*!\brief Evaluation of the expression type type of a subvector operation.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a subvector operation.
// Given the dense or sparse vector type \a VT and the alignment flag \a AF, the nested type
// \a Type corresponds to the resulting return type. In case the given type is neither a
// dense nor a sparse vector type, the resulting data type \a Type is set to \a INVALID_TYPE.
*/
template< typename VT  // Type of the vector operand
        , bool AF >    // Alignment Flag
struct SubvectorExprTrait
{
 private:
   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct DenseResult**************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T >
   struct DenseResult { typedef DenseSubvector<T,AF,IsRowVector<T>::value>  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct SparseResult*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T >
   struct SparseResult { typedef SparseSubvector<T,AF,IsRowVector<T>::value>  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename RemoveReference<VT>::Type  Tmp;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< Or< IsComputation<Tmp>, IsTransExpr<Tmp> >
                      , typename If< Or< IsConst<Tmp>, IsVolatile<Tmp> >
                                   , SubvectorExprTrait< typename RemoveCV<Tmp>::Type, AF >
                                   , Failure
                                   >::Type
                      , typename If< IsDenseVector<Tmp>
                                   , DenseResult<Tmp>
                                   , typename If< IsSparseVector<Tmp>
                                                , SparseResult<Tmp>
                                                , Failure
                                                >::Type
                                   >::Type
                      >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
