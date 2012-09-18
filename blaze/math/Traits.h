//=================================================================================================
/*!
//  \file blaze/math/Traits.h
//  \brief Header file for all expression traits
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

#ifndef _BLAZE_MATH_TRAITS_H_
#define _BLAZE_MATH_TRAITS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/AbsExprTrait.h>
#include <blaze/math/traits/AddExprTrait.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/DMatAbsTrait.h>
#include <blaze/math/traits/DMatDMatAddTrait.h>
#include <blaze/math/traits/DMatDMatMultTrait.h>
#include <blaze/math/traits/DMatDMatSubTrait.h>
#include <blaze/math/traits/DMatDVecMultTrait.h>
#include <blaze/math/traits/DMatScalarDivTrait.h>
#include <blaze/math/traits/DMatScalarMultTrait.h>
#include <blaze/math/traits/DMatSMatAddTrait.h>
#include <blaze/math/traits/DMatSMatMultTrait.h>
#include <blaze/math/traits/DMatSMatSubTrait.h>
#include <blaze/math/traits/DMatSVecMultTrait.h>
#include <blaze/math/traits/DMatTDMatAddTrait.h>
#include <blaze/math/traits/DMatTDMatMultTrait.h>
#include <blaze/math/traits/DMatTDMatSubTrait.h>
#include <blaze/math/traits/DMatTSMatAddTrait.h>
#include <blaze/math/traits/DMatTSMatMultTrait.h>
#include <blaze/math/traits/DMatTSMatSubTrait.h>
#include <blaze/math/traits/DVecAbsTrait.h>
#include <blaze/math/traits/DVecDVecAddTrait.h>
#include <blaze/math/traits/DVecDVecMultTrait.h>
#include <blaze/math/traits/DVecDVecSubTrait.h>
#include <blaze/math/traits/DVecScalarDivTrait.h>
#include <blaze/math/traits/DVecScalarMultTrait.h>
#include <blaze/math/traits/DVecSVecAddTrait.h>
#include <blaze/math/traits/DVecSVecMultTrait.h>
#include <blaze/math/traits/DVecSVecSubTrait.h>
#include <blaze/math/traits/DVecTDVecMultTrait.h>
#include <blaze/math/traits/DVecTSVecMultTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/SMatAbsTrait.h>
#include <blaze/math/traits/SMatDMatAddTrait.h>
#include <blaze/math/traits/SMatDMatMultTrait.h>
#include <blaze/math/traits/SMatDMatSubTrait.h>
#include <blaze/math/traits/SMatDVecMultTrait.h>
#include <blaze/math/traits/SMatScalarDivTrait.h>
#include <blaze/math/traits/SMatScalarMultTrait.h>
#include <blaze/math/traits/SMatSMatAddTrait.h>
#include <blaze/math/traits/SMatSMatMultTrait.h>
#include <blaze/math/traits/SMatSMatSubTrait.h>
#include <blaze/math/traits/SMatSVecMultTrait.h>
#include <blaze/math/traits/SMatTDMatAddTrait.h>
#include <blaze/math/traits/SMatTDMatMultTrait.h>
#include <blaze/math/traits/SMatTDMatSubTrait.h>
#include <blaze/math/traits/SMatTSMatAddTrait.h>
#include <blaze/math/traits/SMatTSMatMultTrait.h>
#include <blaze/math/traits/SMatTSMatSubTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SVecAbsTrait.h>
#include <blaze/math/traits/SVecDVecAddTrait.h>
#include <blaze/math/traits/SVecDVecMultTrait.h>
#include <blaze/math/traits/SVecDVecSubTrait.h>
#include <blaze/math/traits/SVecScalarDivTrait.h>
#include <blaze/math/traits/SVecScalarMultTrait.h>
#include <blaze/math/traits/SVecSVecAddTrait.h>
#include <blaze/math/traits/SVecSVecMultTrait.h>
#include <blaze/math/traits/SVecSVecSubTrait.h>
#include <blaze/math/traits/SVecTDVecMultTrait.h>
#include <blaze/math/traits/SVecTSVecMultTrait.h>
#include <blaze/math/traits/TDMatAbsTrait.h>
#include <blaze/math/traits/TDMatDMatAddTrait.h>
#include <blaze/math/traits/TDMatDMatMultTrait.h>
#include <blaze/math/traits/TDMatDMatSubTrait.h>
#include <blaze/math/traits/TDMatScalarDivTrait.h>
#include <blaze/math/traits/TDMatScalarMultTrait.h>
#include <blaze/math/traits/TDMatSMatAddTrait.h>
#include <blaze/math/traits/TDMatSMatMultTrait.h>
#include <blaze/math/traits/TDMatSMatSubTrait.h>
#include <blaze/math/traits/TDMatDVecMultTrait.h>
#include <blaze/math/traits/TDMatSVecMultTrait.h>
#include <blaze/math/traits/TDMatTDMatAddTrait.h>
#include <blaze/math/traits/TDMatTDMatMultTrait.h>
#include <blaze/math/traits/TDMatTDMatSubTrait.h>
#include <blaze/math/traits/TDMatTSMatAddTrait.h>
#include <blaze/math/traits/TDMatTSMatMultTrait.h>
#include <blaze/math/traits/TDMatTSMatSubTrait.h>
#include <blaze/math/traits/TDVecAbsTrait.h>
#include <blaze/math/traits/TDVecDMatMultTrait.h>
#include <blaze/math/traits/TDVecDVecMultTrait.h>
#include <blaze/math/traits/TDVecScalarDivTrait.h>
#include <blaze/math/traits/TDVecScalarMultTrait.h>
#include <blaze/math/traits/TDVecSMatMultTrait.h>
#include <blaze/math/traits/TDVecSVecMultTrait.h>
#include <blaze/math/traits/TDVecTDMatMultTrait.h>
#include <blaze/math/traits/TDVecTDVecAddTrait.h>
#include <blaze/math/traits/TDVecTDVecMultTrait.h>
#include <blaze/math/traits/TDVecTDVecSubTrait.h>
#include <blaze/math/traits/TDVecTSMatMultTrait.h>
#include <blaze/math/traits/TDVecTSVecAddTrait.h>
#include <blaze/math/traits/TDVecTSVecMultTrait.h>
#include <blaze/math/traits/TDVecTSVecSubTrait.h>
#include <blaze/math/traits/TSMatAbsTrait.h>
#include <blaze/math/traits/TSMatDMatAddTrait.h>
#include <blaze/math/traits/TSMatDMatMultTrait.h>
#include <blaze/math/traits/TSMatDMatSubTrait.h>
#include <blaze/math/traits/TSMatScalarDivTrait.h>
#include <blaze/math/traits/TSMatScalarMultTrait.h>
#include <blaze/math/traits/TSMatSMatAddTrait.h>
#include <blaze/math/traits/TSMatSMatMultTrait.h>
#include <blaze/math/traits/TSMatSMatSubTrait.h>
#include <blaze/math/traits/TSMatDVecMultTrait.h>
#include <blaze/math/traits/TSMatSVecMultTrait.h>
#include <blaze/math/traits/TSMatTDMatAddTrait.h>
#include <blaze/math/traits/TSMatTDMatMultTrait.h>
#include <blaze/math/traits/TSMatTDMatSubTrait.h>
#include <blaze/math/traits/TSMatTSMatAddTrait.h>
#include <blaze/math/traits/TSMatTSMatMultTrait.h>
#include <blaze/math/traits/TSMatTSMatSubTrait.h>
#include <blaze/math/traits/TSVecAbsTrait.h>
#include <blaze/math/traits/TSVecDMatMultTrait.h>
#include <blaze/math/traits/TSVecDVecMultTrait.h>
#include <blaze/math/traits/TSVecScalarDivTrait.h>
#include <blaze/math/traits/TSVecScalarMultTrait.h>
#include <blaze/math/traits/TSVecSMatMultTrait.h>
#include <blaze/math/traits/TSVecSVecMultTrait.h>
#include <blaze/math/traits/TSVecTDMatMultTrait.h>
#include <blaze/math/traits/TSVecTDVecAddTrait.h>
#include <blaze/math/traits/TSVecTDVecMultTrait.h>
#include <blaze/math/traits/TSVecTDVecSubTrait.h>
#include <blaze/math/traits/TSVecTSMatMultTrait.h>
#include <blaze/math/traits/TSVecTSVecAddTrait.h>
#include <blaze/math/traits/TSVecTSVecMultTrait.h>
#include <blaze/math/traits/TSVecTSVecSubTrait.h>

#endif
