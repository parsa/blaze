//=================================================================================================
/*!
//  \file blaze/system/Thresholds.h
//  \brief Header file for the thresholds for matrix/vector and matrix/matrix multiplications
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

#ifndef _BLAZE_SYSTEM_THRESHOLDS_H_
#define _BLAZE_SYSTEM_THRESHOLDS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/StaticAssert.h>




//=================================================================================================
//
//  THRESHOLDS
//
//=================================================================================================

#include <blaze/config/Thresholds.h>


namespace blaze {

//=================================================================================================
//
//  SMP THRESHOLDS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief SMP dense vector assignment threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel assignment of a plain
// dense vector. In case the number of elements of the target vector is larger or equal to this
// threshold, the operation is executed in parallel. If the number of elements is below this
// threshold the operation is executed single-threaded.
*/
const size_t SMP_DVECASSIGN_THRESHOLD = OPENMP_DVECASSIGN_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/dense vector addition threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/dense
// vector addition. In case the number of elements of the target vector is larger or equal to
// this threshold, the operation is executed in parallel. If the number of elements is below
// this threshold the operation is executed single-threaded.
*/
const size_t SMP_DVECDVECADD_THRESHOLD = OPENMP_DVECDVECADD_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/dense vector subtraction threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/dense
// vector subtraction. In case the number of elements of the target vector is larger or equal
// to this threshold, the operation is executed in parallel. If the number of elements is below
// this threshold the operation is executed single-threaded.
*/
const size_t SMP_DVECDVECSUB_THRESHOLD = OPENMP_DVECDVECSUB_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/dense vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/dense
// vector multiplication. In case the number of elements of the target vector is larger or equal
// to this threshold, the operation is executed in parallel. If the number of elements is below
// this threshold the operation is executed single-threaded.
*/
const size_t SMP_DVECDVECMULT_THRESHOLD = OPENMP_DVECDVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/scalar multiplication/division threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/scalar
// multiplication/division. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements
// is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_DVECSCALARMULT_THRESHOLD = OPENMP_DVECSCALARMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/dense vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense matrix/dense
// vector multiplication. In case the number of elements of the target vector is larger or equal
// to this threshold, the operation is executed in parallel. If the number of elements is below
// this threshold the operation is executed single-threaded.
*/
const size_t SMP_DMATDVECMULT_THRESHOLD = OPENMP_DMATDVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major dense matrix/dense vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major dense
// matrix/dense vector multiplication. In case the number of elements of the target vector is
// larger or equal to this threshold, the operation is executed in parallel. If the number of
// elements is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDMATDVECMULT_THRESHOLD = OPENMP_TDMATDVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/row-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/row-major
// dense matrix multiplication. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements is
// below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDVECDMATMULT_THRESHOLD = OPENMP_TDVECDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/column-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/column-major
// dense matrix multiplication. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements is
// below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDVECTDMATMULT_THRESHOLD = OPENMP_TDVECTDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/sparse vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/sparse vector multiplication. In case the number of elements of the target vector is
// larger or equal to this threshold, the operation is executed in parallel. If the number of
// elements is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_DMATSVECMULT_THRESHOLD = OPENMP_DMATSVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major dense matrix/sparse vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major dense
// matrix/sparse vector multiplication. In case the number of elements of the target vector is
// larger or equal to this threshold, the operation is executed in parallel. If the number of
// elements is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDMATSVECMULT_THRESHOLD = OPENMP_TDMATSVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP sparse vector/row-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel sparse vector/row-major
// dense matrix multiplication. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements is
// below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSVECDMATMULT_THRESHOLD = OPENMP_TSVECDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP sparse vector/column-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel sparse vector/column-major
// dense matrix multiplication. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements is
// below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSVECTDMATMULT_THRESHOLD = OPENMP_TSVECTDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major sparse matrix/dense vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major sparse
// matrix/dense vector multiplication. In case the number of elements of the target vector is
// larger or equal to this threshold, the operation is executed in parallel. If the number of
// elements is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_SMATDVECMULT_THRESHOLD = OPENMP_SMATDVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major sparse matrix/dense vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major sparse
// matrix/dense vector multiplication. In case the number of elements of the target vector is
// larger or equal to this threshold, the operation is executed in parallel. If the number of
// elements is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSMATDVECMULT_THRESHOLD = OPENMP_TSMATDVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/row-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/row-major
// sparse matrix multiplication. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements is
// below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDVECSMATMULT_THRESHOLD = OPENMP_TDVECSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/column-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/column-major
// sparse matrix multiplication. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements is
// below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDVECTSMATMULT_THRESHOLD = OPENMP_TDVECTSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major sparse matrix/sparse vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major sparse
// matrix/sparse vector multiplication. In case the number of elements of the target vector is
// larger or equal to this threshold, the operation is executed in parallel. If the number of
// elements is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_SMATSVECMULT_THRESHOLD = OPENMP_SMATSVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major sparse matrix/sparse vector multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major sparse
// matrix/sparse vector multiplication. In case the number of elements of the target vector is
// larger or equal to this threshold, the operation is executed in parallel. If the number of
// elements is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSMATSVECMULT_THRESHOLD = OPENMP_TSMATSVECMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP sparse vector/row-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel sparse vector/row-major
// sparse matrix multiplication. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements is
// below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSVECSMATMULT_THRESHOLD = OPENMP_TSVECSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP sparse vector/column-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel sparse vector/column-major
// sparse matrix multiplication. In case the number of elements of the target vector is larger or
// equal to this threshold, the operation is executed in parallel. If the number of elements is
// below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSVECTSMATMULT_THRESHOLD = OPENMP_TSVECTSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense matrix assignment threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel assignment with a plain
// dense matrix. In case the number of rows/columns of the target matrix is larger or equal to this
// threshold, the operation is executed in parallel. If the number of rows/columns is below this
// threshold the operation is executed single-threaded.
*/
const size_t SMP_DMATASSIGN_THRESHOLD = OPENMP_DMATASSIGN_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/row-major dense matrix addition threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/row-major dense matrix addition. This threshold affects both additions between two
// row-major matrices or two column-major dense matrices. In case the number of rows/columns
// of the target matrix is larger or equal to this threshold, the operation is executed in
// parallel. If the number of rows/columns is below this threshold the operation is executed
// single-threaded.
*/
const size_t SMP_DMATDMATADD_THRESHOLD = OPENMP_DMATDMATADD_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/column-major dense matrix addition threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/column-major dense matrix addition. This threshold affects both additions between a
// row-major matrix and a column-major matrix and a column-major matrix and a row-major matrix.
// In case the number of rows/columns of the target matrix is larger or equal to this threshold,
// the operation is executed in parallel. If the number of rows/columns is below this threshold
// the operation is executed single-threaded.
*/
const size_t SMP_DMATTDMATADD_THRESHOLD = OPENMP_DMATTDMATADD_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/row-major dense matrix subtraction threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/row-major dense matrix subtraction. This threshold affects both subtractions between
// two row-major matrices or two column-major dense matrices. In case the number of rows/columns
// of the target matrix is larger or equal to this threshold, the operation is executed in
// parallel. If the number of rows/columns is below this threshold the operation is executed
// single-threaded.
*/
const size_t SMP_DMATDMATSUB_THRESHOLD = OPENMP_DMATDMATSUB_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/column-major dense matrix subtraction threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/column-major dense matrix subtraction. This threshold affects both subtractions between
// a row-major matrix and a column-major matrix and a column-major matrix and a row-major matrix.
// In case the number of rows/columns of the target matrix is larger or equal to this threshold,
// the operation is executed in parallel. If the number of rows/columns is below this threshold
// the operation is executed single-threaded.
*/
const size_t SMP_DMATTDMATSUB_THRESHOLD = OPENMP_DMATTDMATSUB_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense matrix/scalar multiplication/division threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense matrix/scalar
// multiplication or division can be executed in parallel. In case the number of rows/columns of
// the target matrix is larger or equal to this threshold, the operation is executed in parallel.
// If the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_DMATSCALARMULT_THRESHOLD = OPENMP_DMATSCALARMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/row-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/row-major dense matrix multiplication. In case the number of rows/columns of the target
// matrix is larger or equal to this threshold, the operation is executed in parallel. If the
// number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_DMATDMATMULT_THRESHOLD = OPENMP_DMATDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/column-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/column-major dense matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel.
// If the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_DMATTDMATMULT_THRESHOLD = OPENMP_DMATTDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major dense matrix/row-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major dense
// matrix/row-major dense matrix multiplication. In case the number of rows/columns of the target
// matrix is larger or equal to this threshold, the operation is executed in parallel. If the
// number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDMATDMATMULT_THRESHOLD = OPENMP_TDMATDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major dense matrix/column-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major dense
// matrix/column-major dense matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDMATTDMATMULT_THRESHOLD = OPENMP_TDMATTDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/row-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/row-major sparse matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_DMATSMATMULT_THRESHOLD = OPENMP_DMATSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major dense matrix/column-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major dense
// matrix/column-major sparse matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_DMATTSMATMULT_THRESHOLD = OPENMP_DMATTSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major dense matrix/row-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major dense
// matrix/row-major sparse matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDMATSMATMULT_THRESHOLD = OPENMP_TDMATSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major dense matrix/column-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major dense
// matrix/column-major sparse matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TDMATTSMATMULT_THRESHOLD = OPENMP_TDMATTSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major sparse matrix/row-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major sparse
// matrix/row-major dense matrix multiplication. In case the number of rows/columns of the target
// matrix is larger or equal to this threshold, the operation is executed in parallel. If the
// number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_SMATDMATMULT_THRESHOLD = OPENMP_SMATDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major sparse matrix/column-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major sparse
// matrix/column-major dense matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_SMATTDMATMULT_THRESHOLD = OPENMP_SMATTDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major sparse matrix/row-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major sparse
// matrix/row-major dense matrix multiplication. In case the number of rows/columns of the target
// matrix is larger or equal to this threshold, the operation is executed in parallel. If the
// number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSMATDMATMULT_THRESHOLD = OPENMP_TSMATDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major sparse matrix/column-major dense matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major sparse
// matrix/column-major dense matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel.
// If the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSMATTDMATMULT_THRESHOLD = OPENMP_TSMATTDMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major sparse matrix/row-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major sparse
// matrix/row-major sparse matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_SMATSMATMULT_THRESHOLD = OPENMP_SMATSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP row-major sparse matrix/column-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel row-major sparse
// matrix/column-major sparse matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_SMATTSMATMULT_THRESHOLD = OPENMP_SMATTSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major sparse matrix/row-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major sparse
// matrix/row-major sparse matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSMATSMATMULT_THRESHOLD = OPENMP_TSMATSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP column-major sparse matrix/column-major sparse matrix multiplication threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel column-major sparse
// matrix/column-major sparse matrix multiplication. In case the number of rows/columns of the
// target matrix is larger or equal to this threshold, the operation is executed in parallel. If
// the number of rows/columns is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_TSMATTSMATMULT_THRESHOLD = OPENMP_TSMATTSMATMULT_THRESHOLD;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SMP dense vector/dense vector outer product threshold.
// \ingroup system
//
// This threshold represents the system-specific threshold for a parallel dense vector/dense
// vector outer product. In case the number of rows/columns of the target matrix is larger or
// equal to this threshold, the operation is executed in parallel. If the number of rows/columns
// is below this threshold the operation is executed single-threaded.
*/
const size_t SMP_DVECTDVECMULT_THRESHOLD = OPENMP_DVECTDVECMULT_THRESHOLD;
//*************************************************************************************************

} // namespace blaze




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( blaze::DMATDVECMULT_THRESHOLD   > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDMATDVECMULT_THRESHOLD  > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDVECDMATMULT_THRESHOLD  > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDVECTDMATMULT_THRESHOLD > 0UL );
BLAZE_STATIC_ASSERT( blaze::DMATDMATMULT_THRESHOLD   > 0UL );
BLAZE_STATIC_ASSERT( blaze::DMATTDMATMULT_THRESHOLD  > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDMATDMATMULT_THRESHOLD  > 0UL );
BLAZE_STATIC_ASSERT( blaze::TDMATTDMATMULT_THRESHOLD > 0UL );

BLAZE_STATIC_ASSERT( blaze::OPENMP_DVECASSIGN_THRESHOLD     >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DVECDVECADD_THRESHOLD    >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DVECDVECSUB_THRESHOLD    >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DVECDVECMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DVECSCALARMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATDVECMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDMATDVECMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDVECDMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDVECTDMATMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATSVECMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDMATSVECMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSVECDMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSVECTDMATMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_SMATDVECMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSMATDVECMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDVECSMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDVECTSMATMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_SMATSVECMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSMATSVECMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSVECSMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSVECTSMATMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATASSIGN_THRESHOLD     >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATDMATADD_THRESHOLD    >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATTDMATADD_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATDMATSUB_THRESHOLD    >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATTDMATSUB_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATSCALARMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATDMATMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATTDMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDMATDMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDMATTDMATMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATSMATMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DMATTSMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDMATSMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TDMATTSMATMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_SMATDMATMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_SMATTDMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSMATDMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSMATTDMATMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_SMATSMATMULT_THRESHOLD   >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_SMATTSMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSMATSMATMULT_THRESHOLD  >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_TSMATTSMATMULT_THRESHOLD >= 0UL );
BLAZE_STATIC_ASSERT( blaze::OPENMP_DVECTDVECMULT_THRESHOLD  >= 0UL );

}
/*! \endcond */
//*************************************************************************************************

#endif
