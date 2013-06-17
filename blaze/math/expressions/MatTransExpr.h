//=================================================================================================
/*!
//  \file blaze/math/expressions/MatTransExpr.h
//  \brief Header file for the MatTransExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_MATTRANSEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_MATTRANSEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/TransExpr.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all matrix transposition expression templates.
// \ingroup math
//
// The MatTransExpr class serves as a tag for all expression templates that implement a matrix
// transposition operation. All classes, that represent a matrix transposition operation and
// that are used within the expression template environment of the Blaze library have to derive
// from this class in order to qualify as matrix transposition expression template. Only in case
// a class is derived from the MatTransExpr base class, the IsMatTransExpr type trait recognizes
// the class as valid matrix transposition expression template.
*/
struct MatTransExpr : private TransExpr
{};
//*************************************************************************************************

} // namespace blaze

#endif
