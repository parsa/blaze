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
#include <blaze/math/traits/DMatDMatAddTrait.h>
#include <blaze/math/traits/DMatSMatAddTrait.h>
#include <blaze/math/traits/DMatTDMatAddTrait.h>
#include <blaze/math/traits/DMatTSMatAddTrait.h>
#include <blaze/math/traits/DVecDVecAddTrait.h>
#include <blaze/math/traits/DVecSVecAddTrait.h>
#include <blaze/math/traits/SMatDMatAddTrait.h>
#include <blaze/math/traits/SMatSMatAddTrait.h>
#include <blaze/math/traits/SMatTDMatAddTrait.h>
#include <blaze/math/traits/SMatTSMatAddTrait.h>
#include <blaze/math/traits/SVecDVecAddTrait.h>
#include <blaze/math/traits/SVecSVecAddTrait.h>
#include <blaze/math/traits/TDMatDMatAddTrait.h>
#include <blaze/math/traits/TDMatSMatAddTrait.h>
#include <blaze/math/traits/TDMatTDMatAddTrait.h>
#include <blaze/math/traits/TDMatTSMatAddTrait.h>
#include <blaze/math/traits/TDVecTDVecAddTrait.h>
#include <blaze/math/traits/TDVecTSVecAddTrait.h>
#include <blaze/math/traits/TSMatDMatAddTrait.h>
#include <blaze/math/traits/TSMatSMatAddTrait.h>
#include <blaze/math/traits/TSMatTDMatAddTrait.h>
#include <blaze/math/traits/TSMatTSMatAddTrait.h>
#include <blaze/math/traits/TSVecTDVecAddTrait.h>
#include <blaze/math/traits/TSVecTSVecAddTrait.h>
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
                                                                                       , DMatDMatAddTrait<T1,T2>
                                                                                       , DMatTDMatAddTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TDMatDMatAddTrait<T1,T2>
                                                                                       , TDMatTDMatAddTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , DMatSMatAddTrait<T1,T2>
                                                                                       , DMatTSMatAddTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TDMatSMatAddTrait<T1,T2>
                                                                                       , TDMatTSMatAddTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             >::Type
                                                , typename If< IsDenseMatrix<T2>
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , SMatDMatAddTrait<T1,T2>
                                                                                       , SMatTDMatAddTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TSMatDMatAddTrait<T1,T2>
                                                                                       , TSMatTDMatAddTrait<T1,T2>
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsRowMajorMatrix<T1>
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , SMatSMatAddTrait<T1,T2>
                                                                                       , SMatTSMatAddTrait<T1,T2>
                                                                                       >::Type
                                                                          , typename If< IsRowMajorMatrix<T2>
                                                                                       , TSMatSMatAddTrait<T1,T2>
                                                                                       , TSMatTSMatAddTrait<T1,T2>
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
                                                                                                    , TDVecTDVecAddTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , DVecDVecAddTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TDVecTSVecAddTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , DVecSVecAddTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          >::Type
                                                             , typename If< IsDenseVector<T2>
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TSVecTDVecAddTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , SVecDVecAddTrait<T1,T2>
                                                                                                    >::Type
                                                                                       >::Type
                                                                          , typename If< IsTransposeVector<T1>
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , TSVecTSVecAddTrait<T1,T2>
                                                                                                    , Failure
                                                                                                    >::Type
                                                                                       , typename If< IsTransposeVector<T2>
                                                                                                    , Failure
                                                                                                    , SVecSVecAddTrait<T1,T2>
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
