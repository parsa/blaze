//=================================================================================================
/*!
//  \file blaze/math/Traits.h
//  \brief Header file for all expression traits
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

#ifndef _BLAZE_MATH_TRAITS_H_
#define _BLAZE_MATH_TRAITS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/AddExprTrait.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BinaryMapExprTrait.h>
#include <blaze/math/traits/BinaryMapTrait.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/CrossExprTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/CTransExprTrait.h>
#include <blaze/math/traits/DeclDiagExprTrait.h>
#include <blaze/math/traits/DeclHermExprTrait.h>
#include <blaze/math/traits/DeclLowExprTrait.h>
#include <blaze/math/traits/DeclSymExprTrait.h>
#include <blaze/math/traits/DeclUppExprTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/DMatCTransExprTrait.h>
#include <blaze/math/traits/DMatDeclDiagExprTrait.h>
#include <blaze/math/traits/DMatDeclHermExprTrait.h>
#include <blaze/math/traits/DMatDeclLowExprTrait.h>
#include <blaze/math/traits/DMatDeclSymExprTrait.h>
#include <blaze/math/traits/DMatDeclUppExprTrait.h>
#include <blaze/math/traits/DMatDMatMapExprTrait.h>
#include <blaze/math/traits/DMatEvalExprTrait.h>
#include <blaze/math/traits/DMatMapExprTrait.h>
#include <blaze/math/traits/DMatInvExprTrait.h>
#include <blaze/math/traits/DMatScalarDivExprTrait.h>
#include <blaze/math/traits/DMatScalarMultExprTrait.h>
#include <blaze/math/traits/DMatSerialExprTrait.h>
#include <blaze/math/traits/DMatTDMatMapExprTrait.h>
#include <blaze/math/traits/DMatTransExprTrait.h>
#include <blaze/math/traits/DVecCTransExprTrait.h>
#include <blaze/math/traits/DVecDVecDivExprTrait.h>
#include <blaze/math/traits/DVecDVecMapExprTrait.h>
#include <blaze/math/traits/DVecEvalExprTrait.h>
#include <blaze/math/traits/DVecMapExprTrait.h>
#include <blaze/math/traits/DVecScalarDivExprTrait.h>
#include <blaze/math/traits/DVecScalarMultExprTrait.h>
#include <blaze/math/traits/DVecSerialExprTrait.h>
#include <blaze/math/traits/DVecTransExprTrait.h>
#include <blaze/math/traits/EvalExprTrait.h>
#include <blaze/math/traits/ImagTrait.h>
#include <blaze/math/traits/InvExprTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RealTrait.h>
#include <blaze/math/traits/RowExprTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SchurExprTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SerialExprTrait.h>
#include <blaze/math/traits/SMatCTransExprTrait.h>
#include <blaze/math/traits/SMatDeclDiagExprTrait.h>
#include <blaze/math/traits/SMatDeclHermExprTrait.h>
#include <blaze/math/traits/SMatDeclLowExprTrait.h>
#include <blaze/math/traits/SMatDeclSymExprTrait.h>
#include <blaze/math/traits/SMatDeclUppExprTrait.h>
#include <blaze/math/traits/SMatEvalExprTrait.h>
#include <blaze/math/traits/SMatMapExprTrait.h>
#include <blaze/math/traits/SMatScalarDivExprTrait.h>
#include <blaze/math/traits/SMatScalarMultExprTrait.h>
#include <blaze/math/traits/SMatSerialExprTrait.h>
#include <blaze/math/traits/SMatTransExprTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/traits/SVecCTransExprTrait.h>
#include <blaze/math/traits/SVecDVecDivExprTrait.h>
#include <blaze/math/traits/SVecEvalExprTrait.h>
#include <blaze/math/traits/SVecMapExprTrait.h>
#include <blaze/math/traits/SVecScalarDivExprTrait.h>
#include <blaze/math/traits/SVecScalarMultExprTrait.h>
#include <blaze/math/traits/SVecSerialExprTrait.h>
#include <blaze/math/traits/SVecTransExprTrait.h>
#include <blaze/math/traits/TDMatCTransExprTrait.h>
#include <blaze/math/traits/TDMatDeclDiagExprTrait.h>
#include <blaze/math/traits/TDMatDeclHermExprTrait.h>
#include <blaze/math/traits/TDMatDeclLowExprTrait.h>
#include <blaze/math/traits/TDMatDeclSymExprTrait.h>
#include <blaze/math/traits/TDMatDeclUppExprTrait.h>
#include <blaze/math/traits/TDMatDMatMapExprTrait.h>
#include <blaze/math/traits/TDMatEvalExprTrait.h>
#include <blaze/math/traits/TDMatMapExprTrait.h>
#include <blaze/math/traits/TDMatInvExprTrait.h>
#include <blaze/math/traits/TDMatScalarDivExprTrait.h>
#include <blaze/math/traits/TDMatScalarMultExprTrait.h>
#include <blaze/math/traits/TDMatSerialExprTrait.h>
#include <blaze/math/traits/TDMatTDMatMapExprTrait.h>
#include <blaze/math/traits/TDMatTransExprTrait.h>
#include <blaze/math/traits/TDVecCTransExprTrait.h>
#include <blaze/math/traits/TDVecEvalExprTrait.h>
#include <blaze/math/traits/TDVecMapExprTrait.h>
#include <blaze/math/traits/TDVecScalarDivExprTrait.h>
#include <blaze/math/traits/TDVecScalarMultExprTrait.h>
#include <blaze/math/traits/TDVecSerialExprTrait.h>
#include <blaze/math/traits/TDVecTDVecDivExprTrait.h>
#include <blaze/math/traits/TDVecTDVecMapExprTrait.h>
#include <blaze/math/traits/TDVecTransExprTrait.h>
#include <blaze/math/traits/TransExprTrait.h>
#include <blaze/math/traits/TSMatCTransExprTrait.h>
#include <blaze/math/traits/TSMatDeclDiagExprTrait.h>
#include <blaze/math/traits/TSMatDeclHermExprTrait.h>
#include <blaze/math/traits/TSMatDeclLowExprTrait.h>
#include <blaze/math/traits/TSMatDeclSymExprTrait.h>
#include <blaze/math/traits/TSMatDeclUppExprTrait.h>
#include <blaze/math/traits/TSMatEvalExprTrait.h>
#include <blaze/math/traits/TSMatMapExprTrait.h>
#include <blaze/math/traits/TSMatScalarDivExprTrait.h>
#include <blaze/math/traits/TSMatScalarMultExprTrait.h>
#include <blaze/math/traits/TSMatSerialExprTrait.h>
#include <blaze/math/traits/TSMatTransExprTrait.h>
#include <blaze/math/traits/TSVecCTransExprTrait.h>
#include <blaze/math/traits/TSVecEvalExprTrait.h>
#include <blaze/math/traits/TSVecMapExprTrait.h>
#include <blaze/math/traits/TSVecScalarDivExprTrait.h>
#include <blaze/math/traits/TSVecScalarMultExprTrait.h>
#include <blaze/math/traits/TSVecSerialExprTrait.h>
#include <blaze/math/traits/TSVecTDVecDivExprTrait.h>
#include <blaze/math/traits/TSVecTransExprTrait.h>
#include <blaze/math/traits/UnaryMapExprTrait.h>
#include <blaze/math/traits/UnaryMapTrait.h>

#endif
