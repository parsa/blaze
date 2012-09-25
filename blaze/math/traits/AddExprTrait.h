//=================================================================================================
/*!
//  \file blaze/math/traits/AddExprTrait.h
//  \brief Header file for the AddExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_ADDEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_ADDEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/DMatDMatAddExprTrait.h>
#include <blaze/math/traits/DMatSMatAddExprTrait.h>
#include <blaze/math/traits/DMatTDMatAddExprTrait.h>
#include <blaze/math/traits/DMatTSMatAddExprTrait.h>
#include <blaze/math/traits/DVecDVecAddExprTrait.h>
#include <blaze/math/traits/DVecSVecAddExprTrait.h>
#include <blaze/math/traits/SMatDMatAddExprTrait.h>
#include <blaze/math/traits/SMatSMatAddExprTrait.h>
#include <blaze/math/traits/SMatTDMatAddExprTrait.h>
#include <blaze/math/traits/SMatTSMatAddExprTrait.h>
#include <blaze/math/traits/SVecDVecAddExprTrait.h>
#include <blaze/math/traits/SVecSVecAddExprTrait.h>
#include <blaze/math/traits/TDMatDMatAddExprTrait.h>
#include <blaze/math/traits/TDMatSMatAddExprTrait.h>
#include <blaze/math/traits/TDMatTDMatAddExprTrait.h>
#include <blaze/math/traits/TDMatTSMatAddExprTrait.h>
#include <blaze/math/traits/TDVecTDVecAddExprTrait.h>
#include <blaze/math/traits/TDVecTSVecAddExprTrait.h>
#include <blaze/math/traits/TSMatDMatAddExprTrait.h>
#include <blaze/math/traits/TSMatSMatAddExprTrait.h>
#include <blaze/math/traits/TSMatTDMatAddExprTrait.h>
#include <blaze/math/traits/TSMatTSMatAddExprTrait.h>
#include <blaze/math/traits/TSVecTDVecAddExprTrait.h>
#include <blaze/math/traits/TSVecTSVecAddExprTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
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
/*!\brief Evaluation of the return type of an addition expression.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of an addition expression
// between scalars, vectors, and matrices. Given the two types \a T1 and \a T2, which must
// be either scalar, vector, or matrix types, the nested type \a Type corresponds to the
// resulting return type. In case the types of \a T1 or \a T2 don't fit or if the two types
// cannot be added, the resulting data type \a Type is set to \a INVALID_TYPE.
*/
template< typename T1    // Type of the left-hand side addition operand
        , typename T2 >  // Type of the right-hand side addition operand
struct AddExprTrait
{
 private:
   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { qualified = IsConst<T1>::value || IsVolatile<T1>::value || IsReference<T1>::value ||
                      IsConst<T2>::value || IsVolatile<T2>::value || IsReference<T2>::value };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsMatrix<T1>
                      , typename If< IsMatrix<T2>
                                   , typename If< IsDenseMatrix<T1>
                                                , typename If< IsDenseMatrix<T2>
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , DMatDMatAddExprTrait<T1,T2>
                                                                                       , DMatTDMatAddExprTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TDMatDMatAddExprTrait<T1,T2>
                                                                                       , TDMatTDMatAddExprTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , DMatSMatAddExprTrait<T1,T2>
                                                                                       , DMatTSMatAddExprTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TDMatSMatAddExprTrait<T1,T2>
                                                                                       , TDMatTSMatAddExprTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                , typename If< IsDenseMatrix<T2>
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , SMatDMatAddExprTrait<T1,T2>
                                                                                       , SMatTDMatAddExprTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TSMatDMatAddExprTrait<T1,T2>
                                                                                       , TSMatTDMatAddExprTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , SMatSMatAddExprTrait<T1,T2>
                                                                                       , SMatTSMatAddExprTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TSMatSMatAddExprTrait<T1,T2>
                                                                                       , TSMatTSMatAddExprTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                >::Type
                                   , Failure
                                   >::Type
                      , typename If< IsVector<T1>
                                   , typename If< IsVector<T2>
                                                , typename If< IsDenseVector<T1>
                                                             , typename If< IsDenseVector<T2>
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TDVecTDVecAddExprTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , DVecDVecAddExprTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TDVecTSVecAddExprTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , DVecSVecAddExprTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsDenseVector<T2>
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TSVecTDVecAddExprTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , SVecDVecAddExprTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TSVecTSVecAddExprTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , SVecSVecAddExprTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                , Failure
                                                >::Type
                                   , typename If< IsNumeric<T1>
                                                , typename If< IsNumeric<T2>
                                                             , AddTrait<T1,T2>
                                                             , Failure
                                                             >::Type
                                                , Failure
                                                >::Type
                                   >::Type
                      >::Type  Tmp;

   typedef typename RemoveReference< typename RemoveCV<T1>::Type >::Type  Type1;
   typedef typename RemoveReference< typename RemoveCV<T2>::Type >::Type  Type2;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, AddExprTrait<Type1,Type2>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
