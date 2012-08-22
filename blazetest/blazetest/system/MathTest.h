//=================================================================================================
/*!
//  \file blazemark/system/MathTest.h
//  \brief General settings for the math tests of the blaze test suite
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

#ifndef _BLAZEMARK_SYSTEM_MATHTEST_H_
#define _BLAZEMARK_SYSTEM_MATHTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

//=================================================================================================
//
//  USING DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
using blaze::rowVector;
using blaze::columnVector;
using blaze::rowMajor;
using blaze::columnMajor;

using blaze::CompressedMatrix;
using blaze::CompressedVector;
using blaze::DynamicMatrix;
using blaze::DynamicVector;
using blaze::StaticMatrix;
using blaze::StaticVector;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GENERAL CONFIGURATION
//
//=================================================================================================

#include <blazetest/config/MathTest.h>




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( blaze::IsNumeric<TypeA>::value || blaze::IsVector<TypeA>::value || blaze::IsMatrix<TypeA>::value );
BLAZE_STATIC_ASSERT( blaze::IsNumeric<TypeB>::value || blaze::IsVector<TypeB>::value || blaze::IsMatrix<TypeB>::value );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 2 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION > 2 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION > 2 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_ABS_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_ABS_OPERATION > 2 ) );

}
/*! \endcond */
//*************************************************************************************************

} // namespace mathtest

} // namespace blazetest

#endif
