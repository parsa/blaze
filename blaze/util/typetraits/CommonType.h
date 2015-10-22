//=================================================================================================
/*!
//  \file blaze/util/typetraits/CommonType.h
//  \brief Header file for the CommonType type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_COMMONTYPE_H_
#define _BLAZE_UTIL_TYPETRAITS_COMMONTYPE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/common_type.hpp>
#include <blaze/util/NullType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Deduction of a type common to several types.
// \ingroup type_traits
//
// The CommonType type trait deduces a type that is common to up to three types \a T1 to
// \a T3. This may for instance be the resulting type of an arithmetic operation, such as an
// addition or a subtraction. Note that cv and reference qualifiers are generally ignored.

   \code
   blaze::CommonType<short,int>::Type                         // Results in 'int'
   blaze::CommonType<const double,int&>::Type                 // Results in 'double'
   blaze::CommonType<char&, volatile int, const float>::Type  // Results in 'float'
   \endcode
*/
template< typename T1, typename T2, typename T3 = NullType >
struct CommonType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename boost::common_type<T1,T2,T3>::type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Partial specialization of the CommonType type trait for two given types.
template< typename T1, typename T2 >
struct CommonType<T1,T2,NullType>
{
 public:
   //**********************************************************************************************
   typedef typename boost::common_type<T1,T2>::type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
