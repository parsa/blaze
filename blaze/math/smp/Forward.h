//=================================================================================================
/*!
//  \file blaze/math/smp/Forward.h
//  \brief Header file for all SMP forward declarations
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

#ifndef _BLAZE_MATH_SMP_FORWARD_H_
#define _BLAZE_MATH_SMP_FORWARD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Forward.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/util/EnableIf.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> >;

template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpAddAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> >;

template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpSubAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> >;

template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpMultAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> >;

template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpDivAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> >;


template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsSparseVector_v<VT1> >;

template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpAddAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsSparseVector_v<VT1> >;

template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpSubAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsSparseVector_v<VT1> >;

template< typename VT1, bool TF1, typename VT2, bool TF2 >
auto smpMultAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsSparseVector_v<VT1> >;


template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpAddAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpSubAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpSchurAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >;


template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsSparseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpAddAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsSparseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpSubAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsSparseMatrix_v<MT1> >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
auto smpSchurAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsSparseMatrix_v<MT1> >;

} // namespace blaze

#endif
