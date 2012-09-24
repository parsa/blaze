//=================================================================================================
/*!
//  \file blaze/math/traits/MultExprTrait.h
//  \brief Header file for the MultExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_MULTEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_MULTEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/DMatDMatMultTrait.h>
#include <blaze/math/traits/DMatDVecMultTrait.h>
#include <blaze/math/traits/DMatScalarMultTrait.h>
#include <blaze/math/traits/DMatSMatMultTrait.h>
#include <blaze/math/traits/DMatSVecMultTrait.h>
#include <blaze/math/traits/DMatTDMatMultTrait.h>
#include <blaze/math/traits/DMatTSMatMultTrait.h>
#include <blaze/math/traits/DVecDVecMultTrait.h>
#include <blaze/math/traits/DVecScalarMultTrait.h>
#include <blaze/math/traits/DVecSVecMultTrait.h>
#include <blaze/math/traits/DVecTDVecMultTrait.h>
#include <blaze/math/traits/DVecTSVecMultTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SMatDMatMultTrait.h>
#include <blaze/math/traits/SMatDVecMultTrait.h>
#include <blaze/math/traits/SMatScalarMultTrait.h>
#include <blaze/math/traits/SMatSMatMultTrait.h>
#include <blaze/math/traits/SMatSVecMultTrait.h>
#include <blaze/math/traits/SMatTDMatMultTrait.h>
#include <blaze/math/traits/SMatTSMatMultTrait.h>
#include <blaze/math/traits/SVecDVecMultTrait.h>
#include <blaze/math/traits/SVecScalarMultTrait.h>
#include <blaze/math/traits/SVecSVecMultTrait.h>
#include <blaze/math/traits/SVecTDVecMultTrait.h>
#include <blaze/math/traits/SVecTSVecMultTrait.h>
#include <blaze/math/traits/TDMatDMatMultTrait.h>
#include <blaze/math/traits/TDMatDVecMultTrait.h>
#include <blaze/math/traits/TDMatScalarMultTrait.h>
#include <blaze/math/traits/TDMatSMatMultTrait.h>
#include <blaze/math/traits/TDMatSVecMultTrait.h>
#include <blaze/math/traits/TDMatTDMatMultTrait.h>
#include <blaze/math/traits/TDMatTSMatMultTrait.h>
#include <blaze/math/traits/TDVecDMatMultTrait.h>
#include <blaze/math/traits/TDVecDVecMultTrait.h>
#include <blaze/math/traits/TDVecScalarMultTrait.h>
#include <blaze/math/traits/TDVecSMatMultTrait.h>
#include <blaze/math/traits/TDVecSVecMultTrait.h>
#include <blaze/math/traits/TDVecTDMatMultTrait.h>
#include <blaze/math/traits/TDVecTDVecMultTrait.h>
#include <blaze/math/traits/TDVecTSMatMultTrait.h>
#include <blaze/math/traits/TDVecTSVecMultTrait.h>
#include <blaze/math/traits/TSMatDMatMultTrait.h>
#include <blaze/math/traits/TSMatDVecMultTrait.h>
#include <blaze/math/traits/TSMatScalarMultTrait.h>
#include <blaze/math/traits/TSMatSMatMultTrait.h>
#include <blaze/math/traits/TSMatSVecMultTrait.h>
#include <blaze/math/traits/TSMatTDMatMultTrait.h>
#include <blaze/math/traits/TSMatTSMatMultTrait.h>
#include <blaze/math/traits/TSVecDMatMultTrait.h>
#include <blaze/math/traits/TSVecDVecMultTrait.h>
#include <blaze/math/traits/TSVecScalarMultTrait.h>
#include <blaze/math/traits/TSVecSMatMultTrait.h>
#include <blaze/math/traits/TSVecSVecMultTrait.h>
#include <blaze/math/traits/TSVecTDMatMultTrait.h>
#include <blaze/math/traits/TSVecTDVecMultTrait.h>
#include <blaze/math/traits/TSVecTSMatMultTrait.h>
#include <blaze/math/traits/TSVecTSVecMultTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/IfNot.h>
#include <blaze/util/SelectType.h>
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
/*!\brief Evaluation of the resulting expression type of a multiplication.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a multiplication expression
// between scalars, vectors, and matrices. Given the two types \a T1 and \a T2, which must be
// either scalar, vector, or matrix types, the nested type \a Type corresponds to the resulting
// return type. In case \a T1 or \a T2 don't fit or if the two types cannot be multiplied, the
// resulting data type \a Type is set to \a INVALID_TYPE.
*/
template< typename T1    // Type of the left-hand side multiplication operand
        , typename T2 >  // Type of the right-hand side multiplication operand
struct MultExprTrait
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
                      , typename If< IsDenseMatrix<T1>
                                   , typename If< IsRowMajorMatrix<T1>
                                                , typename If< IsMatrix<T2>
                                                             , typename If< IsDenseMatrix<T2>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , DMatDMatMultTrait<T1,T2>
                                                                                       , DMatTDMatMultTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , DMatSMatMultTrait<T1,T2>
                                                                                       , DMatTSMatMultTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsVector<T2>
                                                                          , typename If< IsDenseVector<T2>
                                                                                       , typename IfNot< IsTransposeVector<T2>
                                                                                                       , DMatDVecMultTrait<T1,T2>
                                                                                                       , Failure
                                                                                                       >::Type
                                                                                       , typename IfNot< IsTransposeVector<T2>
                                                                                                       , DMatSVecMultTrait<T1,T2>
                                                                                                       , Failure
                                                                                                       >::Type
                                                                                       >::Type
                                                                          , typename If< IsNumeric<T2>
                                                                                       , DMatScalarMultTrait<T1,T2>
                                                                                       , Failure
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                , typename If< IsMatrix<T2>
                                                             , typename If< IsDenseMatrix<T2>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TDMatDMatMultTrait<T1,T2>
                                                                                       , TDMatTDMatMultTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TDMatSMatMultTrait<T1,T2>
                                                                                       , TDMatTSMatMultTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsVector<T2>
                                                                          , typename If< IsDenseVector<T2>
                                                                                       , typename IfNot< IsTransposeVector<T2>
                                                                                                       , TDMatDVecMultTrait<T1,T2>
                                                                                                       , Failure
                                                                                                       >::Type
                                                                                       , typename IfNot< IsTransposeVector<T2>
                                                                                                       , TDMatSVecMultTrait<T1,T2>
                                                                                                       , Failure
                                                                                                       >::Type
                                                                                       >::Type
                                                                          , typename If< IsNumeric<T2>
                                                                                       , TDMatScalarMultTrait<T1,T2>
                                                                                       , Failure
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                >::Type
                                   , typename If< IsRowMajorMatrix<T1>
                                                , typename If< IsMatrix<T2>
                                                             , typename If< IsDenseMatrix<T2>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , SMatDMatMultTrait<T1,T2>
                                                                                       , SMatTDMatMultTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , SMatSMatMultTrait<T1,T2>
                                                                                       , SMatTSMatMultTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsVector<T2>
                                                                          , typename If< IsDenseVector<T2>
                                                                                       , typename IfNot< IsTransposeVector<T2>
                                                                                                       , SMatDVecMultTrait<T1,T2>
                                                                                                       , Failure
                                                                                                       >::Type
                                                                                       , typename IfNot< IsTransposeVector<T2>
                                                                                                       , SMatSVecMultTrait<T1,T2>
                                                                                                       , Failure
                                                                                                       >::Type
                                                                                       >::Type
                                                                          , typename If< IsNumeric<T2>
                                                                                       , SMatScalarMultTrait<T1,T2>
                                                                                       , Failure
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                , typename If< IsMatrix<T2>
                                                             , typename If< IsDenseMatrix<T2>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TSMatDMatMultTrait<T1,T2>
                                                                                       , TSMatTDMatMultTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TSMatSMatMultTrait<T1,T2>
                                                                                       , TSMatTSMatMultTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsVector<T2>
                                                                          , typename If< IsDenseVector<T2>
                                                                                       , typename IfNot< IsTransposeVector<T2>
                                                                                                       , TSMatDVecMultTrait<T1,T2>
                                                                                                       , Failure
                                                                                                       >::Type
                                                                                       , typename IfNot< IsTransposeVector<T2>
                                                                                                       , TSMatSVecMultTrait<T1,T2>
                                                                                                       , Failure
                                                                                                       >::Type
                                                                                       >::Type
                                                                          , typename If< IsNumeric<T2>
                                                                                       , TSMatScalarMultTrait<T1,T2>
                                                                                       , Failure
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                >::Type
                                   >::Type
                      , typename If< IsVector<T1>
                                   , typename If< IsDenseVector<T1>
                                                , typename If< IsTransposeVector<T1>
                                                             , typename If< IsMatrix<T2>
                                                                          , typename If< IsDenseMatrix<T2>
                                                                                       , typename If< IsRowMajorMatrix<T2>
                                                                                                    , TDVecDMatMultTrait<T1,T2>
                                                                                                    , TDVecTDMatMultTrait<T1,T2>
                                                                                                    >::Type
                                                                                       , typename If< IsRowMajorMatrix<T2>
                                                                                                    , TDVecSMatMultTrait<T1,T2>
                                                                                                    , TDVecTSMatMultTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsVector<T2>
                                                                                       , typename If< IsDenseVector<T2>
                                                                                                    , typename If< IsTransposeVector<T2>
                                                                                                                 , TDVecTDVecMultTrait<T1,T2>
                                                                                                                 , TDVecDVecMultTrait<T1,T2>
                                                                                                                 >::Type
                                                                                                    , typename If< IsTransposeVector<T2>
                                                                                                                 , TDVecTSVecMultTrait<T1,T2>
                                                                                                                 , TDVecSVecMultTrait<T1,T2>
                                                                                                                 >::Type
                                                                                                    >::Type
                                                                                       , typename If< IsNumeric<T2>
                                                                                                    , TDVecScalarMultTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsVector<T2>
                                                                          , typename If< IsDenseVector<T2>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , DVecTDVecMultTrait<T1,T2>
                                                                                                    , DVecDVecMultTrait<T1,T2>
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , DVecTSVecMultTrait<T1,T2>
                                                                                                    , DVecSVecMultTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsNumeric<T2>
                                                                                       , DVecScalarMultTrait<T1,T2>
                                                                                       , Failure
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                , typename If< IsTransposeVector<T1>
                                                             , typename If< IsMatrix<T2>
                                                                          , typename If< IsDenseMatrix<T2>
                                                                                       , typename If< IsRowMajorMatrix<T2>
                                                                                                    , TSVecDMatMultTrait<T1,T2>
                                                                                                    , TSVecTDMatMultTrait<T1,T2>
                                                                                                    >::Type
                                                                                       , typename If< IsRowMajorMatrix<T2>
                                                                                                    , TSVecSMatMultTrait<T1,T2>
                                                                                                    , TSVecTSMatMultTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsVector<T2>
                                                                                       , typename If< IsDenseVector<T2>
                                                                                                    , typename If< IsTransposeVector<T2>
                                                                                                                 , TSVecTDVecMultTrait<T1,T2>
                                                                                                                 , TSVecDVecMultTrait<T1,T2>
                                                                                                                 >::Type
                                                                                                    , typename If< IsTransposeVector<T2>
                                                                                                                 , TSVecTSVecMultTrait<T1,T2>
                                                                                                                 , TSVecSVecMultTrait<T1,T2>
                                                                                                                 >::Type
                                                                                                    >::Type
                                                                                       , typename If< IsNumeric<T2>
                                                                                                    , TSVecScalarMultTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsVector<T2>
                                                                          , typename If< IsDenseVector<T2>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , SVecTDVecMultTrait<T1,T2>
                                                                                                    , SVecDVecMultTrait<T1,T2>
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , SVecTSVecMultTrait<T1,T2>
                                                                                                    , SVecSVecMultTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsNumeric<T2>
                                                                                       , SVecScalarMultTrait<T1,T2>
                                                                                       , Failure
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                >::Type
                                   , typename If< IsNumeric<T1>
                                                , typename If< IsMatrix<T2>
                                                    , typename If< IsDenseMatrix<T2>
                                                                 , typename If< IsRowMajorMatrix<T2>
                                                                              , DMatScalarMultTrait<T2,T1>
                                                                              , TDMatScalarMultTrait<T2,T1>
                                                                              >::Type
                                                                 , typename If< IsRowMajorMatrix<T2>
                                                                              , SMatScalarMultTrait<T2,T1>
                                                                              , TSMatScalarMultTrait<T2,T1>
                                                                              >::Type
                                                                 >::Type
                                                    , typename If< IsVector<T2>
                                                                 , typename If< IsDenseVector<T2>
                                                                              , typename If< IsTransposeVector<T2>
                                                                                           , TDVecScalarMultTrait<T2,T1>
                                                                                           , DVecScalarMultTrait<T2,T1>
                                                                                           >::Type
                                                                              , typename If< IsTransposeVector<T2>
                                                                                           , TSVecScalarMultTrait<T2,T1>
                                                                                           , SVecScalarMultTrait<T2,T1>
                                                                                           >::Type
                                                                              >::Type
                                                                 , typename If< IsNumeric<T2>
                                                                              , MultTrait<T1,T2>
                                                                              , Failure
                                                                              >::Type
                                                                 >::Type
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
   typedef typename SelectType< qualified, MultExprTrait<Type1,Type2>, Tmp >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
