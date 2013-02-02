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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Assert.h>
#include <blaze/util/logging/FunctionTrace.h>


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




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Vector global functions */
//@{
template< typename VT, bool TF >
inline size_t size( const Vector<VT,TF>& v );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline void assign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline void addAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline void subAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline void multAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current size/dimension of the vector.
// \ingroup vector
//
// \param v The given vector.
// \return The size of the vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
inline size_t size( Vector<VT,TF>& v )
{
   return (~v).size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default assignment of a vector to another vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline void assign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;
   
   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );
   (~lhs).assign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be added.
// \return void
//
// This function implements the default addition assignment of a vector to a vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline void addAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;
   
   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );
   (~lhs).addAssign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be subtracted.
// \return void
//
// This function implements the default subtraction assignment of a vector to a vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline void subAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;
   
   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );
   (~lhs).subAssign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the multiplication assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be multiplied.
// \return void
//
// This function implements the default multiplication assignment of a vector to a vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline void multAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;
   
   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );
   (~lhs).multAssign( ~rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
