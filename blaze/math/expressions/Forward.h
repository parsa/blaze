//=================================================================================================
/*!
//  \file blaze/math/expressions/Forward.h
//  \brief Header file for all forward declarations for expression class templates
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
template< typename, bool > class DMatConjExpr;
template< typename, typename, bool > class DMatDMatAddExpr;
template< typename, typename > class DMatDMatMultExpr;
template< typename, typename, bool > class DMatDMatSubExpr;
template< typename, typename > class DMatDVecMultExpr;
template< typename, bool > class DMatEvalExpr;
template< typename, bool > class DMatImagExpr;
template< typename, bool > class DMatInvExpr;
template< typename, bool > class DMatRealExpr;
template< typename, typename, bool > class DMatScalarDivExpr;
template< typename, typename, bool > class DMatScalarMultExpr;
template< typename, bool > class DMatSerialExpr;
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
template< typename, bool > class DVecConjExpr;
template< typename, typename, bool > class DVecDVecAddExpr;
template< typename, typename > class DVecDVecCrossExpr;
template< typename, typename, bool > class DVecDVecMultExpr;
template< typename, typename, bool > class DVecDVecSubExpr;
template< typename, bool > class DVecEvalExpr;
template< typename, bool > class DVecImagExpr;
template< typename, bool > class DVecRealExpr;
template< typename, typename, bool > class DVecScalarDivExpr;
template< typename, typename, bool > class DVecScalarMultExpr;
template< typename, bool > class DVecSerialExpr;
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
template< typename, bool > class SMatConjExpr;
template< typename, typename > class SMatDMatMultExpr;
template< typename, typename, bool > class SMatDMatSubExpr;
template< typename, typename > class SMatDVecMultExpr;
template< typename, bool > class SMatEvalExpr;
template< typename, bool > class SMatImagExpr;
template< typename, bool > class SMatRealExpr;
template< typename, typename, bool > class SMatScalarDivExpr;
template< typename, typename, bool > class SMatScalarMultExpr;
template< typename, bool > class SMatSerialExpr;
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
template< typename, bool > class SVecConjExpr;
template< typename, typename > class SVecDVecCrossExpr;
template< typename, typename, bool > class SVecDVecMultExpr;
template< typename, typename, bool > class SVecDVecSubExpr;
template< typename, bool > class SVecEvalExpr;
template< typename, bool > class SVecImagExpr;
template< typename, bool > class SVecRealExpr;
template< typename, typename, bool > class SVecScalarDivExpr;
template< typename, typename, bool > class SVecScalarMultExpr;
template< typename, bool > class SVecSerialExpr;
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

template< typename VT, bool TF >
inline const DVecTransExpr<VT,!TF> trans( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
inline const SVecTransExpr<VT,!TF> trans( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
inline const DMatTransExpr<MT,!SO> trans( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
inline const SMatTransExpr<MT,!SO> trans( const SparseMatrix<MT,SO>& );

template< typename VT, bool TF >
inline const DVecSerialExpr<VT,TF> serial( const DenseVector<VT,TF>& );

template< typename VT, bool TF >
inline const SVecSerialExpr<VT,TF> serial( const SparseVector<VT,TF>& );

template< typename MT, bool SO >
inline const DMatSerialExpr<MT,SO> serial( const DenseMatrix<MT,SO>& );

template< typename MT, bool SO >
inline const SMatSerialExpr<MT,SO> serial( const SparseMatrix<MT,SO>& );

} // namespace blaze

#endif
