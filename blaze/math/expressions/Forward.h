//=================================================================================================
/*!
//  \file blaze/math/expressions/Forward.h
//  \brief Header file for all forward declarations for expression class templates
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

#ifndef _BLAZE_MATH_EXPRESSIONS_FORWARD_H_
#define _BLAZE_MATH_EXPRESSIONS_FORWARD_H_


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< typename, bool > struct DenseMatrix;
template< typename, bool > struct DenseVector;
template< typename, bool > class DMatAbsExpr;
template< typename, typename, bool > class DMatDMatAddExpr;
template< typename, typename > class DMatDMatMultExpr;
template< typename, typename, bool > class DMatDMatSubExpr;
template< typename, typename > class DMatDVecMultExpr;
template< typename, bool > class DMatEvalExpr;
template< typename, typename, bool > class DMatScalarDivExpr;
template< typename, typename, bool > class DMatScalarMultExpr;
template< typename, typename, bool > class DMatSMatAddExpr;
template< typename, typename > class DMatSMatMultExpr;
template< typename, typename, bool > class DMatSMatSubExpr;
template< typename, typename > class DMatSVecMultExpr;
template< typename, typename > class DMatTDMatAddExpr;
template< typename, typename > class DMatTDMatMultExpr;
template< typename, typename > class DMatTDMatSubExpr;
template< typename, bool > class DMatTransExpr;
template< typename, bool > class DMatTransposer;
template< typename, typename > class DMatTSMatAddExpr;
template< typename, typename > class DMatTSMatMultExpr;
template< typename, typename > class DMatTSMatSubExpr;
template< typename, bool > class DVecAbsExpr;
template< typename, typename, bool > class DVecDVecAddExpr;
template< typename, typename > class DVecDVecCrossExpr;
template< typename, typename, bool > class DVecDVecMultExpr;
template< typename, typename, bool > class DVecDVecSubExpr;
template< typename, bool > class DVecEvalExpr;
template< typename, typename, bool > class DVecScalarDivExpr;
template< typename, typename, bool > class DVecScalarMultExpr;
template< typename, typename, bool > class DVecSVecAddExpr;
template< typename, typename > class DVecSVecCrossExpr;
template< typename, typename, bool > class DVecSVecMultExpr;
template< typename, typename, bool > class DVecSVecSubExpr;
template< typename, typename > class DVecTDVecMultExpr;
template< typename, bool > class DVecTransExpr;
template< typename, bool > class DVecTransposer;
template< typename, typename > class DVecTSVecMultExpr;
template< typename, bool > struct Matrix;
template< typename, bool > class SMatAbsExpr;
template< typename, typename > class SMatDMatMultExpr;
template< typename, typename, bool > class SMatDMatSubExpr;
template< typename, typename > class SMatDVecMultExpr;
template< typename, bool > class SMatEvalExpr;
template< typename, typename, bool > class SMatScalarDivExpr;
template< typename, typename, bool > class SMatScalarMultExpr;
template< typename, typename > class SMatSMatAddExpr;
template< typename, typename > class SMatSMatMultExpr;
template< typename, typename > class SMatSMatSubExpr;
template< typename, typename > class SMatSVecMultExpr;
template< typename, typename > class SMatTDMatMultExpr;
template< typename, typename > class SMatTDMatSubExpr;
template< typename, bool > class SMatTransExpr;
template< typename, bool > class SMatTransposer;
template< typename, typename > class SMatTSMatAddExpr;
template< typename, typename > class SMatTSMatMultExpr;
template< typename, typename > class SMatTSMatSubExpr;
template< typename, bool > struct SparseMatrix;
template< typename, bool > struct SparseVector;
template< typename, bool > class SVecAbsExpr;
template< typename, typename > class SVecDVecCrossExpr;
template< typename, typename, bool > class SVecDVecMultExpr;
template< typename, typename, bool > class SVecDVecSubExpr;
template< typename, bool > class SVecEvalExpr;
template< typename, typename, bool > class SVecScalarDivExpr;
template< typename, typename, bool > class SVecScalarMultExpr;
template< typename, typename, bool > class SVecSVecAddExpr;
template< typename, typename > class SVecSVecCrossExpr;
template< typename, typename, bool > class SVecSVecMultExpr;
template< typename, typename, bool > class SVecSVecSubExpr;
template< typename, typename > class SVecTDVecMultExpr;
template< typename, bool > class SVecTransExpr;
template< typename, bool > class SVecTransposer;
template< typename, typename > class SVecTSVecMultExpr;
template< typename, typename > class TDMatDMatMultExpr;
template< typename, typename > class TDMatDVecMultExpr;
template< typename, typename > class TDMatSMatAddExpr;
template< typename, typename > class TDMatSMatMultExpr;
template< typename, typename > class TDMatSMatSubExpr;
template< typename, typename > class TDMatSVecMultExpr;
template< typename, typename > class TDMatTDMatMultExpr;
template< typename, typename > class TDMatTSMatMultExpr;
template< typename, typename > class TDVecDMatMultExpr;
template< typename, typename > class TDVecSMatMultExpr;
template< typename, typename > class TDVecTDMatMultExpr;
template< typename, typename > class TDVecTSMatMultExpr;
template< typename, typename > class TSMatDMatMultExpr;
template< typename, typename > class TSMatDMatSubExpr;
template< typename, typename > class TSMatDVecMultExpr;
template< typename, typename > class TSMatSMatMultExpr;
template< typename, typename > class TSMatSMatSubExpr;
template< typename, typename > class TSMatSVecMultExpr;
template< typename, typename > class TSMatTDMatMultExpr;
template< typename, typename > class TSMatTSMatAddExpr;
template< typename, typename > class TSMatTSMatMultExpr;
template< typename, typename > class TSMatTSMatSubExpr;
template< typename, typename > class TSVecDMatMultExpr;
template< typename, typename > class TSVecSMatMultExpr;
template< typename, typename > class TSVecTDMatMultExpr;
template< typename, typename > class TSVecTSMatMultExpr;
template< typename, bool > struct Vector;

} // namespace blaze

#endif
