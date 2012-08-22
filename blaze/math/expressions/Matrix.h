//=================================================================================================
/*!
//  \file blaze/math/expressions/Matrix.h
//  \brief Header file for the Matrix base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_MATRIX_H_
#define _BLAZE_MATH_EXPRESSIONS_MATRIX_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup matrix Matrices
// \ingroup math
*/
/*!\brief Base class for matrices.
// \ingroup matrix
//
// The Matrix class is a base class for all dense and sparse matrix classes within the Blaze
// library. It provides an abstraction from the actual type of the matrix, but enables a
// conversion back to this type via the 'Curiously Recurring Template Pattern' (CRTP).
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
struct Matrix
{
   //**Type definitions****************************************************************************
   typedef MT  MatrixType;  //!< Type of the matrix.
   //**********************************************************************************************

   //**Non-const conversion operator***************************************************************
   /*!\brief Conversion operator for non-constant matrices.
   //
   // \return Reference of the actual type of the matrix.
   */
   inline MatrixType& operator~() {
      return *static_cast<MatrixType*>( this );
   }
   //**********************************************************************************************

   //**Const conversion operator*******************************************************************
   /*!\brief Conversion operator for constant matrices.
   //
   // \return Constant reference of the actual type of the matrix.
   */
   inline const MatrixType& operator~() const {
      return *static_cast<const MatrixType*>( this );
   }
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
