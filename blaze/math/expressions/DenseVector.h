//=================================================================================================
/*!
//  \file blaze/math/expressions/DenseVector.h
//  \brief Header file for the DenseVector base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DENSEVECTOR_H_
#define _BLAZE_MATH_EXPRESSIONS_DENSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Vector.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dense_vector Dense Vectors
// \ingroup vector
*/
/*!\defgroup dense_vector_expression Expressions
// \ingroup dense_vector
*/
/*!\brief Base class for N-dimensional dense vectors.
// \ingroup dense_vector
//
// The DenseVector class is a base class for all arbitrarily sized (N-dimensional) dense
// vectors. It provides an abstraction from the actual type of the dense vector, but enables
// a conversion back to this type via the Vector base class.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
struct DenseVector : public Vector<VT,TF>
{};
//*************************************************************************************************

} // namespace blaze

#endif
