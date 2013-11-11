//=================================================================================================
/*!
//  \file blaze/util/SizeTrait.h
//  \brief Header file for the size trait
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

#ifndef _BLAZE_UTIL_SIZETRAIT_H_
#define _BLAZE_UTIL_SIZETRAIT_H_


namespace blaze {

//=================================================================================================
//
//  SIZETRAITHELPER CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary class for the SizeTrait class template.
// \ingroup util
*/
template< typename T1, typename T2, bool >
struct SizeTraitHelper
{
   //**Type definitions****************************************************************************
   typedef T1  Large;
   typedef T2  Small;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the auxiliary SizeTraitHelper class template.
// \ingroup util
*/
template< typename T1, typename T2 >
struct SizeTraitHelper<T1,T2,true>
{
   //**Type definitions****************************************************************************
   typedef T2  Large;
   typedef T1  Small;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZETRAIT CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the SizeTrait class.
// \ingroup util
//
// The SizeTrait class template evaluates the larger of the two given data types by use of the
// sizeof operator. SizeTrait defines the data types \a Large for the larger of the two given
// data types and \a Small for the smaller data type. In case both data types have the same
// size, the first given data type \a T1 is chosen as the large and \a T2 as the small data
// type.
*/
template< typename T1, typename T2 >
struct SizeTrait
{
 private:
   //**Type definitions****************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Instantiation of the SizeTraitHelper class template.
   typedef SizeTraitHelper<T1,T2,( sizeof(T2) > sizeof(T1) )>  Helper;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename Helper::Large  Large;  //!< The larger of the two given data types.
   typedef typename Helper::Small  Small;  //!< The smaller of the two given data types.
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
