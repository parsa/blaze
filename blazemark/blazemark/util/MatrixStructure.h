//=================================================================================================
/*!
//  \file blazemark/util/MatrixStructure.h
//  \brief Header file for the matrix structure flags
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

#ifndef _BLAZEMARK_UTIL_MATRIXSTRUCTURE_H_
#define _BLAZEMARK_UTIL_MATRIXSTRUCTURE_H_


namespace blazemark {

//=================================================================================================
//
//  MATRIX STRUCTURE FLAGS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Matrix structure flags.
//
// The MatrixStructure enumeration represents all possible structures of (sparse) matrices.
*/
enum MatrixStructure {
   band   = 1,  //!< Flag for banded matrices.
                /*!< The \a band flag indicates a banded matrix of a certain bandwidth. The
                     following example demonstrates a typical 3-banded matrix resulting from
                     a physical setting:
                     \f[\left(\begin{array}{*{5}{c}}
                     2  & -1 & 0  & 0  & 0  \\
                     -1 & 2  & -1 & 0  & 0  \\
                     0  & -1 & 2  & -1 & 0  \\
                     0  & 0  & -1 & 2  & -1 \\
                     0  & 0  & 0  & -1 & 2  \\
                     \end{array}\right)\f]. */
   random = 2   //!< Flag for random matrices.
                /*!< The \a random flag indicates a matrix with randomly determined non-zero
                     entries. The following example demonstrates a random matrix with 2
                     non-zero entries per row:
                     \f[\left(\begin{array}{*{5}{c}}
                     0 & 3 & 0 & 0 & 2 \\
                     1 & 0 & 5 & 0 & 0 \\
                     0 & 1 & 7 & 0 & 0 \\
                     0 & 0 & 1 & 0 & 4 \\
                     0 & 8 & 0 & 1 & 0 \\
                     \end{array}\right)\f]. */
};
//*************************************************************************************************

} // namespace blazemark

#endif
