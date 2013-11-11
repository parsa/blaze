//=================================================================================================
/*!
//  \file blaze/util/SelectType.h
//  \brief Header file for the SelectType class template
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

#ifndef _BLAZE_UTIL_SELECTTYPE_H_
#define _BLAZE_UTIL_SELECTTYPE_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type selection.
// \ingroup util
//
// The SelectType class template selects one of the two given types \a T1 and \a T2 depending
// on the \a Select template argument. In case the \a Select compile time constant expression
// evaluates to \a true, the member type definition \a Type is set to \a T1. In case \a Select
// evaluates to \a false, \a Type is set to \a T2.
*/
template< bool Select    // Compile time selection
        , typename T1    // Type to be selected if Select=true
        , typename T2 >  // Type to be selected if Select=false
struct SelectType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef T1  Type;  //!< The selected type.
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the SwitchType class template.
// \ingroup util
//
// This specialization of the SelectType template is selected in case the \a Select compile time
// constant expression evaluates to \a false. The member type definition is set to the second
// given type \a T2.
*/
template< typename T1    // Type not to be selected
        , typename T2 >  // Type to be selected
struct SelectType<false,T1,T2>
{
 public:
   //**********************************************************************************************
   typedef T2  Type;  //!< The selected type.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
