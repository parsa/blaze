//=================================================================================================
/*!
//  \file blaze/math/expressions/MatMatAddExpr.h
//  \brief Header file for the MatMatAddExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_MATMATADDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_MATMATADDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/AddExpr.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all matrix/matrix addition expression templates.
// \ingroup math
//
// The MatMatAddExpr class serves as a tag for all expression templates that implement a
// matrix/matrix addition. All classes, that represent a matrix addition and that are used
// within the expression template environment of the Blaze library have to derive from this
// class in order to qualify as matrix addition expression template. Only in case a class is
// derived from the MatMatAddExpr base class, the IsMatMatAddExpr type trait recognizes the
// class as valid matrix addition expression template.
*/
struct MatMatAddExpr : private AddExpr
{};
//*************************************************************************************************

} // namespace blaze

#endif
