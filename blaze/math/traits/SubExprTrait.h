//=================================================================================================
/*!
//  \file blaze/math/traits/SubExprTrait.h
//  \brief Header file for the SubExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_SUBEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_SUBEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/MathTrait.h>
#include <blaze/math/traits/DMatDMatSubTrait.h>
#include <blaze/math/traits/DMatSMatSubTrait.h>
#include <blaze/math/traits/DMatTDMatSubTrait.h>
#include <blaze/math/traits/DMatTSMatSubTrait.h>
#include <blaze/math/traits/DVecDVecSubTrait.h>
#include <blaze/math/traits/DVecSVecSubTrait.h>
#include <blaze/math/traits/SMatDMatSubTrait.h>
#include <blaze/math/traits/SMatSMatSubTrait.h>
#include <blaze/math/traits/SMatTDMatSubTrait.h>
#include <blaze/math/traits/SMatTSMatSubTrait.h>
#include <blaze/math/traits/SVecDVecSubTrait.h>
#include <blaze/math/traits/SVecSVecSubTrait.h>
#include <blaze/math/traits/TDMatDMatSubTrait.h>
#include <blaze/math/traits/TDMatSMatSubTrait.h>
#include <blaze/math/traits/TDMatTDMatSubTrait.h>
#include <blaze/math/traits/TDMatTSMatSubTrait.h>
#include <blaze/math/traits/TDVecTDVecSubTrait.h>
#include <blaze/math/traits/TDVecTSVecSubTrait.h>
#include <blaze/math/traits/TSMatDMatSubTrait.h>
#include <blaze/math/traits/TSMatSMatSubTrait.h>
#include <blaze/math/traits/TSMatTDMatSubTrait.h>
#include <blaze/math/traits/TSMatTSMatSubTrait.h>
#include <blaze/math/traits/TSVecTDVecSubTrait.h>
#include <blaze/math/traits/TSVecTSVecSubTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the return type of a subtraction expression.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the return type of a subtraction expression
// between scalar, vectors, and matrices. Given the two types \a T1 and \a T2, which must be
// either scalar, vector, or matrix types, the nested type \a Type corresponds to the resulting
// return type. In case \a T1 or \a T2 don't fit or if the two types cannot be subtracted, the
// resulting data type \a Type is set to \a INVALID_TYPE.
*/
template< typename T1    // Type of the left-hand side subtraction operand
        , typename T2 >  // Type of the right-hand side subtraction operand
struct SubExprTrait
{
 private:
   //**struct ScalarSub****************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename ST1, typename ST2 >
   struct ScalarSub { typedef typename MathTrait<ST1,ST2>::SubType  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct Failure******************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< IsMatrix<T1>
                      , typename If< IsMatrix<T2>
                                   , typename If< IsDenseMatrix<T1>
                                                , typename If< IsDenseMatrix<T2>
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , DMatDMatSubTrait<T1,T2>
                                                                                       , DMatTDMatSubTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TDMatDMatSubTrait<T1,T2>
                                                                                       , TDMatTDMatSubTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , DMatSMatSubTrait<T1,T2>
                                                                                       , DMatTSMatSubTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TDMatSMatSubTrait<T1,T2>
                                                                                       , TDMatTSMatSubTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                , typename If< IsDenseMatrix<T2>
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , SMatDMatSubTrait<T1,T2>
                                                                                       , SMatTDMatSubTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TSMatDMatSubTrait<T1,T2>
                                                                                       , TSMatTDMatSubTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , SMatSMatSubTrait<T1,T2>
                                                                                       , SMatTSMatSubTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TSMatSMatSubTrait<T1,T2>
                                                                                       , TSMatTSMatSubTrait<T1,T2>
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
                                                                                                    , TDVecTDVecSubTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , DVecDVecSubTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TDVecTSVecSubTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , DVecSVecSubTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsDenseVector<T2>
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TSVecTDVecSubTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , SVecDVecSubTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TSVecTSVecSubTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , SVecSVecSubTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                , Failure
                                                >::Type
                                   , typename If< IsNumeric<T1>
                                                , typename If< IsNumeric<T2>
                                                             , ScalarSub<T1,T2>
                                                             , Failure
                                                             >::Type
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
