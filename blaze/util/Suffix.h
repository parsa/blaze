//=================================================================================================
/*!
//  \file blaze/util/Suffix.h
//  \brief Header file for compile time constraints
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

#ifndef _BLAZE_UTIL_SUFFIX_H_
#define _BLAZE_UTIL_SUFFIX_H_


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Helper macro for macro concatenation.
// \ingroup util
//
// The following code was borrowed from the Boost C++ framework (www.boost.org). This piece of
// macro magic joins the two arguments together, even when one of the arguments is itself a
// macro (see 16.3.1 in C++ standard).  The key is that macro expansion of macro arguments does
// not occur in BLAZE_DO_JOIN2 but does in BLAZE_DO_JOIN.
*/
#define BLAZE_JOIN( X, Y ) BLAZE_DO_JOIN( X, Y )
#define BLAZE_DO_JOIN( X, Y ) BLAZE_DO_JOIN2(X,Y)
#define BLAZE_DO_JOIN2( X, Y ) X##Y
/*! \endcond */
//*************************************************************************************************

#endif
