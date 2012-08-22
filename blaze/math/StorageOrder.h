//=================================================================================================
/*!
//  \file blaze/math/StorageOrder.h
//  \brief Header file for the matrix storage order types
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

#ifndef _BLAZE_MATH_STORAGEORDER_H_
#define _BLAZE_MATH_STORAGEORDER_H_


namespace blaze {

//=================================================================================================
//
//  MATRIX STORAGE ORDER TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Storage order flag for row-major matrices.
//
// Via this flag it is possible to specify the storage order of matrices as row-major. For
// instance, given the following matrix

                          \f[\left(\begin{array}{*{3}{c}}
                          1 & 2 & 3 \\
                          4 & 5 & 6 \\
                          \end{array}\right)\f]\n

// in case of row-major order the elements are stored in the order

                          \f[\left(\begin{array}{*{6}{c}}
                          1 & 2 & 3 & 4 & 5 & 6. \\
                          \end{array}\right)\f]

// The following example demonstrates the setup of this \f$ 2 \times 3 \f$ matrix:

   \code
   using blaze::rowMajor;
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> A( 1, 2, 3, 4, 5, 6 );
   \endcode
*/
const bool rowMajor = false;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Storage order flag for column-major matrices.
//
// Via this flag it is possible to specify the storage order of matrices as column-major. For
// instance, given the following matrix

                          \f[\left(\begin{array}{*{3}{c}}
                          1 & 2 & 3 \\
                          4 & 5 & 6 \\
                          \end{array}\right)\f]\n

// in case of column-major order the elements are stored in the order

                          \f[\left(\begin{array}{*{6}{c}}
                          1 & 4 & 2 & 5 & 3 & 6. \\
                          \end{array}\right)\f]

// The following example demonstrates the setup of this \f$ 2 \times 3 \f$ matrix:

   \code
   using blaze::columnMajor;
   blaze::StaticMatrix<int,2UL,3UL,columnMajor> A( 1, 2, 3, 4, 5, 6 );
   \endcode
*/
const bool columnMajor = true;
//*************************************************************************************************

} // namespace blaze

#endif
