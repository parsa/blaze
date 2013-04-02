//=================================================================================================
/*!
//  \file blaze/util/AlignmentTrait.h
//  \brief Header file for the alignment trait
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

#ifndef _BLAZE_UTIL_ALIGNMENTTRAIT_H_
#define _BLAZE_UTIL_ALIGNMENTTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/alignment_of.hpp>
#include <blaze/system/Vectorization.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/typetraits/IsVectorizable.h>


namespace blaze {

//=================================================================================================
//
//  SIZETRAIT CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the required alignment of the given data type.
// \ingroup util
//
// The AlignmentTrait class template evaluates the required alignment for the given data type.
// For instance, for fundamental data types that can be vectorized via SSE or AVX instructions,
// the proper alignment is 16 or 32 bytes, respectively. For all other data types, a multiple
// of the alignment chosen by the compiler is returned. The evaluated alignment can be queried
// via the nested \a value member.

   \code
   AlignmentTrait<unsigned int>::value  // Evaluates to 16 if SSE is available, a multiple of
                                        // the alignment chosen by the compiler otherwise.
   AlignmentTrait<double>::value        // Evaluates to 32 if AVX is available, to 16 if only
                                        // SSE is available, and a multiple of the alignment
                                        // chosen by the compiler otherwise.
   \endcode
*/
template< typename T >
struct AlignmentTrait
{
 public:
   //**Member enumerations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
#if BLAZE_MIC_MODE
   enum { value = ( IsVectorizable<T>::value )?( 64UL ):( boost::alignment_of<T>::value ) };
#elif BLAZE_SSE2_MODE
   enum { value = ( IsVectorizable<T>::value )?( 16UL ):( boost::alignment_of<T>::value ) };
#else
   enum { value = boost::alignment_of<T>::value };
#endif
   /*! \endcond */
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST   ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( T );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignmentTrait for 'float'.
// \ingroup util
*/
template<>
struct AlignmentTrait<float>
{
 public:
   //**Member enumerations*************************************************************************
#if BLAZE_MIC_MODE
   enum { value = 64UL };
#elif BLAZE_AVX_MODE
   enum { value = 32UL };
#elif BLAZE_SSE_MODE
   enum { value = 16UL };
#else
   enum { value = boost::alignment_of<float>::value };
#endif
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of AlignmentTrait for 'double'.
// \ingroup util
*/
template<>
struct AlignmentTrait<double>
{
 public:
   //**Member enumerations*************************************************************************
#if BLAZE_MIC_MODE
   enum { value = 64UL };
#elif BLAZE_AVX_MODE
   enum { value = 32UL };
#elif BLAZE_SSE_MODE
   enum { value = 16UL };
#else
   enum { value = boost::alignment_of<double>::value };
#endif
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
