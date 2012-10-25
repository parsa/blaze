//=================================================================================================
/*!
//  \file blaze/math/intrinsics/IntrinsicTrait.h
//  \brief Header file for the intrinsic trait
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

#ifndef _BLAZE_MATH_INTRINSICS_INTRINSICTRAIT_H_
#define _BLAZE_MATH_INTRINSICS_INTRINSICTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_cv.hpp>
#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/SSE.h>
#include <blaze/util/AlignmentTrait.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Volatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS INTRINSICTRAITHELPER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IntrinsicTrait class template.
// \ingroup intrinsics
//
// This helper structure provides the mapping between the size of an integral data type and the
// according intrinsic data type.
*/
template< size_t N >
struct IntrinsicTraitHelper;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitHelper class template for 1-byte integral data types.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitHelper<1UL>
{
   typedef sse_int8_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 16 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = 0 };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitHelper class template for 2-byte integral data types.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitHelper<2UL>
{
   typedef sse_int16_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 8 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = BLAZE_SSE2_MODE };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitHelper class template for 4-byte integral data types.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
struct IntrinsicTraitHelper<4UL>
{
   typedef sse_int32_t  Type;
   enum { size           = 16,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1 };
};
#else
template<>
struct IntrinsicTraitHelper<4UL>
{
   typedef sse_int32_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 4 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = BLAZE_SSE4_MODE };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitHelper class template for 8-byte integral data types.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitHelper<8UL>
{
   typedef sse_int64_t  Type;
   enum { size           = 8,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1 };
};
#else
template<>
struct IntrinsicTraitHelper<8UL>
{
   typedef sse_int64_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 2 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = 0 };
};
#endif
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS INTRINSICTRAITBASE
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Base template for the IntrinsicTraitBase class.
// \ingroup intrinsics
*/
template< typename T >
struct IntrinsicTraitBase
{
   typedef T  Type;
   enum { size           = 1,
          alignment      = AlignmentTrait<T>::value,
          addition       = 0,
          subtraction    = 0,
          multiplication = 0 };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'short'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<short>
{
 private:
   typedef IntrinsicTraitHelper<sizeof(short)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentTrait<short>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'unsigned short'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<unsigned short>
{
 private:
   typedef IntrinsicTraitHelper<sizeof(unsigned short)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentTrait<unsigned short>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'int'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<int>
{
 private:
   typedef IntrinsicTraitHelper<sizeof(int)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentTrait<int>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'unsigned int'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<unsigned int>
{
 private:
   typedef IntrinsicTraitHelper<sizeof(unsigned int)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentTrait<unsigned int>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'long'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<long>
{
 private:
   typedef IntrinsicTraitHelper<sizeof(long)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentTrait<long>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'unsigned long'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<unsigned long>
{
 private:
   typedef IntrinsicTraitHelper<sizeof(unsigned long)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentTrait<unsigned long>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'float'.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
struct IntrinsicTraitBase<float>
{
   typedef sse_float_t  Type;
   enum { size           = ( 64UL / sizeof(float) ),
          alignment      = AlignmentTrait<float>::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1 };
};
#elif BLAZE_AVX_MODE
template<>
struct IntrinsicTraitBase<float>
{
   typedef sse_float_t  Type;
   enum { size           = ( 32UL / sizeof(float) ),
          alignment      = AlignmentTrait<float>::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1 };
};
#else
template<>
struct IntrinsicTraitBase<float>
{
   typedef sse_float_t  Type;
   enum { size           = ( BLAZE_SSE_MODE )?( 16UL / sizeof(float) ):( 1 ),
          alignment      = AlignmentTrait<float>::value,
          addition       = BLAZE_SSE_MODE,
          subtraction    = BLAZE_SSE_MODE,
          multiplication = BLAZE_SSE_MODE };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'double'.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitBase<double>
{
   typedef sse_double_t  Type;
   enum { size           = ( 64UL / sizeof(double) ),
          alignment      = AlignmentTrait<double>::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1 };
};
#elif BLAZE_AVX_MODE
template<>
struct IntrinsicTraitBase<double>
{
   typedef sse_double_t  Type;
   enum { size           = ( 32UL / sizeof(double) ),
          alignment      = AlignmentTrait<double>::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1 };
};
#else
template<>
struct IntrinsicTraitBase<double>
{
   typedef sse_double_t  Type;
   enum { size           = ( BLAZE_SSE_MODE )?( 16UL / sizeof(double) ):( 1 ),
          alignment      = AlignmentTrait<double>::value,
          addition       = BLAZE_SSE_MODE,
          subtraction    = BLAZE_SSE_MODE,
          multiplication = BLAZE_SSE_MODE };
};
#endif
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS INTRINSICTRAIT
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Intrinsic characteristics of data types.
// \ingroup intrinsics
//
// The IntrinsicTrait class template provides the intrinsic characteristics of a specific data
// type:
//
//  - The nested data type \a Type corresponds to the according packed, intrinsic data type. In
//    case the data type doesn't have an intrinsic representation, \a Type corresonds to the given
//    data type itself.
//  - The \a size value corresponds to the number of values of the given data type that are packed
//    together in one intrinsic vector type. In case the data type cannot be vectorized, \a size
//    is set to 1.
//  - If the data type can be involved in vectorized additions, the \a addition value is set to 1.
//    Otherwise, \a addition is set to 0.
//  - In case the data type supports vectorized subtractions, the \a subtraction value is set to 1.
//    Else it is set to 0.
//  - If the data type supports vectorized multiplications, the \a multiplication value is set to
//    1. If it cannot be used in multiplications, it is set to 0.
*/
template< typename T >
class IntrinsicTrait : public IntrinsicTraitBase< typename boost::remove_cv<T>::type >
{
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST   ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( T );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
