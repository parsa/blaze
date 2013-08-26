//=================================================================================================
/*!
//  \file blaze/math/traits/SubvectorExprTrait.h
//  \brief Header file for the SubvectorExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_SUBVECTOREXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_SUBVECTOREXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsTransExpr.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/math/views/Forward.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/StaticAssert.h>
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
// Via this type trait it is possible to evaluate the return type of a subvector operation. Given
// the dense or sparse vector type \a VT, the nested type \a Type corresponds to the resulting
// return type. In case the given type is neither a dense nor a sparse vector type, the resulting
// data type \a Type is set to \a INVALID_TYPE.
*/
template< typename VT >  // Type of the vector operand
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
   struct DenseResult { typedef DenseSubvector<T,IsTransposeVector<T>::value>  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct SparseResult*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename T >
   struct SparseResult { typedef SparseSubvector<T,IsTransposeVector<T>::value>  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename RemoveReference< typename RemoveCV<VT>::Type >::Type  Tmp;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< Or< IsComputation<VT>, IsTransExpr<VT> >
                      , SubvectorExprTrait<Tmp>
                      , typename If< IsDenseVector<VT>
                                   , DenseResult<VT>
                                   , typename If< IsSparseVector<VT>
                                                , SparseResult<VT>
                                                , Failure
                                                >::Type
                                   >::Type
                      >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_STATIC_ASSERT( !( Or< IsComputation<VT>, IsTransExpr<VT> >::value ) ||
                        IsConst<VT>::value || IsVolatile<VT>::value || IsReference<VT>::value );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
