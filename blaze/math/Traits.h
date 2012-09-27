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
#include <blaze/math/traits/CMathTrait.h>
#include <blaze/math/traits/CrossExprTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/DMatAbsExprTrait.h>
#include <blaze/math/traits/DMatDMatAddExprTrait.h>
#include <blaze/math/traits/DMatDMatMultExprTrait.h>
#include <blaze/math/traits/DMatDMatSubExprTrait.h>
#include <blaze/math/traits/DMatDVecMultExprTrait.h>
#include <blaze/math/traits/DMatScalarDivExprTrait.h>
#include <blaze/math/traits/DMatScalarMultExprTrait.h>
#include <blaze/math/traits/DMatSMatAddExprTrait.h>
#include <blaze/math/traits/DMatSMatMultExprTrait.h>
#include <blaze/math/traits/DMatSMatSubExprTrait.h>
#include <blaze/math/traits/DMatSVecMultExprTrait.h>
#include <blaze/math/traits/DMatTDMatAddExprTrait.h>
#include <blaze/math/traits/DMatTDMatMultExprTrait.h>
#include <blaze/math/traits/DMatTDMatSubExprTrait.h>
#include <blaze/math/traits/DMatTSMatAddExprTrait.h>
#include <blaze/math/traits/DMatTSMatMultExprTrait.h>
#include <blaze/math/traits/DMatTSMatSubExprTrait.h>
#include <blaze/math/traits/DVecAbsExprTrait.h>
#include <blaze/math/traits/DVecDVecAddExprTrait.h>
#include <blaze/math/traits/DVecDVecCrossExprTrait.h>
#include <blaze/math/traits/DVecDVecMultExprTrait.h>
#include <blaze/math/traits/DVecDVecSubExprTrait.h>
#include <blaze/math/traits/DVecScalarDivExprTrait.h>
#include <blaze/math/traits/DVecScalarMultExprTrait.h>
#include <blaze/math/traits/DVecSVecAddExprTrait.h>
#include <blaze/math/traits/DVecSVecCrossExprTrait.h>
#include <blaze/math/traits/DVecSVecMultExprTrait.h>
#include <blaze/math/traits/DVecSVecSubExprTrait.h>
#include <blaze/math/traits/DVecTDVecMultExprTrait.h>
#include <blaze/math/traits/DVecTSVecMultExprTrait.h>
#include <blaze/math/traits/MathTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SMatAbsExprTrait.h>
#include <blaze/math/traits/SMatDMatAddExprTrait.h>
#include <blaze/math/traits/SMatDMatMultExprTrait.h>
#include <blaze/math/traits/SMatDMatSubExprTrait.h>
#include <blaze/math/traits/SMatDVecMultExprTrait.h>
#include <blaze/math/traits/SMatScalarDivExprTrait.h>
#include <blaze/math/traits/SMatScalarMultExprTrait.h>
#include <blaze/math/traits/SMatSMatAddExprTrait.h>
#include <blaze/math/traits/SMatSMatMultExprTrait.h>
#include <blaze/math/traits/SMatSMatSubExprTrait.h>
#include <blaze/math/traits/SMatSVecMultExprTrait.h>
#include <blaze/math/traits/SMatTDMatAddExprTrait.h>
#include <blaze/math/traits/SMatTDMatMultExprTrait.h>
#include <blaze/math/traits/SMatTDMatSubExprTrait.h>
#include <blaze/math/traits/SMatTSMatAddExprTrait.h>
#include <blaze/math/traits/SMatTSMatMultExprTrait.h>
#include <blaze/math/traits/SMatTSMatSubExprTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SVecAbsExprTrait.h>
#include <blaze/math/traits/SVecDVecAddExprTrait.h>
#include <blaze/math/traits/SVecDVecCrossExprTrait.h>
#include <blaze/math/traits/SVecDVecMultExprTrait.h>
#include <blaze/math/traits/SVecDVecSubExprTrait.h>
#include <blaze/math/traits/SVecScalarDivExprTrait.h>
#include <blaze/math/traits/SVecScalarMultExprTrait.h>
#include <blaze/math/traits/SVecSVecAddExprTrait.h>
#include <blaze/math/traits/SVecSVecCrossExprTrait.h>
#include <blaze/math/traits/SVecSVecMultExprTrait.h>
#include <blaze/math/traits/SVecSVecSubExprTrait.h>
#include <blaze/math/traits/SVecTDVecMultExprTrait.h>
#include <blaze/math/traits/SVecTSVecMultExprTrait.h>
#include <blaze/math/traits/TDMatAbsExprTrait.h>
#include <blaze/math/traits/TDMatDMatAddExprTrait.h>
#include <blaze/math/traits/TDMatDMatMultExprTrait.h>
#include <blaze/math/traits/TDMatDMatSubExprTrait.h>
#include <blaze/math/traits/TDMatScalarDivExprTrait.h>
#include <blaze/math/traits/TDMatScalarMultExprTrait.h>
#include <blaze/math/traits/TDMatSMatAddExprTrait.h>
#include <blaze/math/traits/TDMatSMatMultExprTrait.h>
#include <blaze/math/traits/TDMatSMatSubExprTrait.h>
#include <blaze/math/traits/TDMatDVecMultExprTrait.h>
#include <blaze/math/traits/TDMatSVecMultExprTrait.h>
#include <blaze/math/traits/TDMatTDMatAddExprTrait.h>
#include <blaze/math/traits/TDMatTDMatMultExprTrait.h>
#include <blaze/math/traits/TDMatTDMatSubExprTrait.h>
#include <blaze/math/traits/TDMatTSMatAddExprTrait.h>
#include <blaze/math/traits/TDMatTSMatMultExprTrait.h>
#include <blaze/math/traits/TDMatTSMatSubExprTrait.h>
#include <blaze/math/traits/TDVecAbsExprTrait.h>
#include <blaze/math/traits/TDVecDMatMultExprTrait.h>
#include <blaze/math/traits/TDVecDVecMultExprTrait.h>
#include <blaze/math/traits/TDVecScalarDivExprTrait.h>
#include <blaze/math/traits/TDVecScalarMultExprTrait.h>
#include <blaze/math/traits/TDVecSMatMultExprTrait.h>
#include <blaze/math/traits/TDVecSVecMultExprTrait.h>
#include <blaze/math/traits/TDVecTDMatMultExprTrait.h>
#include <blaze/math/traits/TDVecTDVecAddExprTrait.h>
#include <blaze/math/traits/TDVecTDVecMultExprTrait.h>
#include <blaze/math/traits/TDVecTDVecSubExprTrait.h>
#include <blaze/math/traits/TDVecTSMatMultExprTrait.h>
#include <blaze/math/traits/TDVecTSVecAddExprTrait.h>
#include <blaze/math/traits/TDVecTSVecMultExprTrait.h>
#include <blaze/math/traits/TDVecTSVecSubExprTrait.h>
#include <blaze/math/traits/TSMatAbsExprTrait.h>
#include <blaze/math/traits/TSMatDMatAddExprTrait.h>
#include <blaze/math/traits/TSMatDMatMultExprTrait.h>
#include <blaze/math/traits/TSMatDMatSubExprTrait.h>
#include <blaze/math/traits/TSMatScalarDivExprTrait.h>
#include <blaze/math/traits/TSMatScalarMultExprTrait.h>
#include <blaze/math/traits/TSMatSMatAddExprTrait.h>
#include <blaze/math/traits/TSMatSMatMultExprTrait.h>
#include <blaze/math/traits/TSMatSMatSubExprTrait.h>
#include <blaze/math/traits/TSMatDVecMultExprTrait.h>
#include <blaze/math/traits/TSMatSVecMultExprTrait.h>
#include <blaze/math/traits/TSMatTDMatAddExprTrait.h>
#include <blaze/math/traits/TSMatTDMatMultExprTrait.h>
#include <blaze/math/traits/TSMatTDMatSubExprTrait.h>
#include <blaze/math/traits/TSMatTSMatAddExprTrait.h>
#include <blaze/math/traits/TSMatTSMatMultExprTrait.h>
#include <blaze/math/traits/TSMatTSMatSubExprTrait.h>
#include <blaze/math/traits/TSVecAbsExprTrait.h>
#include <blaze/math/traits/TSVecDMatMultExprTrait.h>
#include <blaze/math/traits/TSVecDVecMultExprTrait.h>
#include <blaze/math/traits/TSVecScalarDivExprTrait.h>
#include <blaze/math/traits/TSVecScalarMultExprTrait.h>
#include <blaze/math/traits/TSVecSMatMultExprTrait.h>
#include <blaze/math/traits/TSVecSVecMultExprTrait.h>
#include <blaze/math/traits/TSVecTDMatMultExprTrait.h>
#include <blaze/math/traits/TSVecTDVecAddExprTrait.h>
#include <blaze/math/traits/TSVecTDVecMultExprTrait.h>
#include <blaze/math/traits/TSVecTDVecSubExprTrait.h>
#include <blaze/math/traits/TSVecTSMatMultExprTrait.h>
#include <blaze/math/traits/TSVecTSVecAddExprTrait.h>
#include <blaze/math/traits/TSVecTSVecMultExprTrait.h>
#include <blaze/math/traits/TSVecTSVecSubExprTrait.h>

#endif
