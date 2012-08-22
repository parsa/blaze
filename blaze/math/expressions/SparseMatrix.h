//=================================================================================================
/*!
//  \file blaze/math/expressions/SparseMatrix.h
//  \brief Header file for the SparseMatrix base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SPARSEMATRIX_H_
#define _BLAZE_MATH_EXPRESSIONS_SPARSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Matrix.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup sparse_matrix Sparse Matrices
// \ingroup matrix
*/
/*!\defgroup sparse_matrix_expression Expressions
// \ingroup sparse_matrix
*/
/*!\brief Base class for sparse matrices.
// \ingroup sparse_matrix
//
// The SparseMatrix class is a base class for all sparse matrix classes. It provides an
// abstraction from the actual type of the sparse matrix, but enables a conversion back
// to this type via the Matrix base class.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
struct SparseMatrix : public Matrix<MT,SO>
{};
//*************************************************************************************************

} // namespace blaze

#endif
