//=================================================================================================
/*!
//  \file blaze/math/traits/BinaryMapExprTrait.h
//  \brief Header file for the BinaryMapExprTrait class template
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#include <utility>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/RemoveReference.h>


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

   //**struct Result*******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Result { using Type = decltype( map( std::declval<T1>()
                                             , std::declval<T2>()
                                             , std::declval<OP>() ) ); };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_t< ( IsVector_v< RemoveReference_t<T1> > && IsVector_v< RemoveReference_t<T2> > ) ||
                               ( IsMatrix_v< RemoveReference_t<T1> > && IsMatrix_v< RemoveReference_t<T2> > )
                             , Result
                             , Failure >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the BinaryMapExprTrait class template.
// \ingroup math_traits
//
// The BinaryMapExprTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the BinaryMapExprTrait class template. For instance, given the data types \a T1 and
// and \a T2 and the custom operation type \a OP the following two type definitions are identical:

   \code
   using Type1 = typename BinaryMapExprTrait<T1,T2,OP>::Type;
   using Type2 = BinaryMapExprTrait_t<T1,T2,OP>;
   \endcode
*/
template< typename T1    // Type of the left-hand side map operand
        , typename T2    // Type of the right-hand side map operand
        , typename OP >  // Type of the custom operation
using BinaryMapExprTrait_t = typename BinaryMapExprTrait<T1,T2,OP>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
