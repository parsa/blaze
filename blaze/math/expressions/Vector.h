//=================================================================================================
/*!
//  \file blaze/math/expressions/Vector.h
//  \brief Header file for the Vector CRTP base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_VECTOR_H_
#define _BLAZE_MATH_EXPRESSIONS_VECTOR_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup vector Vectors
// \ingroup math
*/
/*!\brief Base class for N-dimensional vectors.
// \ingroup vector
//
// The Vector class is a base class for all arbitrarily sized (N-dimensional) dense and sparse
// vector classes within the Blaze library. It provides an abstraction from the actual type of
// the vector, but enables a conversion back to this type via the 'Curiously Recurring Template
// Pattern' (CRTP).
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
struct Vector
{
   //**Type definitions****************************************************************************
   typedef VT  VectorType;  //!< Type of the vector.
   //**********************************************************************************************

   //**Non-const conversion operator***************************************************************
   /*!\brief Conversion operator for non-constant vectors.
   //
   // \return Reference of the actual type of the vector.
   */
   inline VectorType& operator~() {
      return *static_cast<VectorType*>( this );
   }
   //**********************************************************************************************

   //**Const conversion operators******************************************************************
   /*!\brief Conversion operator for constant vectors.
   //
   // \return Const reference of the actual type of the vector.
   */
   inline const VectorType& operator~() const {
      return *static_cast<const VectorType*>( this );
   }
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
