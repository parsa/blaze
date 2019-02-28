//=================================================================================================
/*!
//  \file blaze/util/DecltypeAuto.h
//  \brief Header file for the decltype(auto) workaround
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_UTIL_DECLTYPEAUTO_H_
#define _BLAZE_UTIL_DECLTYPEAUTO_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Compiler.h>


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compiler specific patch for decltype(auto) variable definitions.
// \ingroup util
//
// The BLAZE_DECLTYPE_AUTO is a patch for the Intel C++ compiler that sometimes emits an internal
// compiler error when encountering a variable definition of the following form:

   \code
   decltype(auto) var( expr );
   \endcode

// In order to circumvent this compilation error, the BLAZE_DECLTYPE_AUTO macro can be used:

   \code
   BLAZE_DECLTYPE_AUTO( var, expr );
   \endcode
*/
#if BLAZE_INTEL_COMPILER
#  define BLAZE_DECLTYPE_AUTO(VAR,EXPR) decltype( EXPR ) VAR( EXPR )
#else
#  define BLAZE_DECLTYPE_AUTO(VAR,EXPR) decltype(auto) VAR( EXPR )
#endif
/*! \endcond */
//*************************************************************************************************

#endif
