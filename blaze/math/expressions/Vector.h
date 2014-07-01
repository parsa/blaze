//=================================================================================================
/*!
//  \file blaze/math/expressions/Vector.h
//  \brief Header file for the Vector CRTP base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_VECTOR_H_
#define _BLAZE_MATH_EXPRESSIONS_VECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Assert.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/typetraits/IsSame.h>


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
inline typename VT::Iterator begin( Vector<VT,TF>& v );

template< typename VT, bool TF >
inline typename VT::ConstIterator begin( const Vector<VT,TF>& v );

template< typename VT, bool TF >
inline typename VT::ConstIterator cbegin( const Vector<VT,TF>& v );

template< typename VT, bool TF >
inline typename VT::Iterator end( Vector<VT,TF>& v );

template< typename VT, bool TF >
inline typename VT::ConstIterator end( const Vector<VT,TF>& v );

template< typename VT, bool TF >
inline typename VT::ConstIterator cend( const Vector<VT,TF>& v );

template< typename VT, bool TF >
inline size_t size( const Vector<VT,TF>& v );

template< typename VT, bool TF >
inline size_t capacity( const Vector<VT,TF>& v );

template< typename VT, bool TF >
inline size_t nonZeros( const Vector<VT,TF>& v );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline void assign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline void addAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline void subAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline void multAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs );

template< typename VT1, bool TF1, typename VT2, bool TF2 >
inline bool isSame( const Vector<VT1,TF1>& a, const Vector<VT2,TF2>& b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given vector.
//
// \param v The given dense or sparse vector.
// \return Iterator to the first element of the given vector.
*/
template< typename VT, bool TF >
inline typename VT::Iterator begin( Vector<VT,TF>& v )
{
   return (~v).begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given vector.
//
// \param v The given dense or sparse vector.
// \return Iterator to the first element of the given vector.
*/
template< typename VT, bool TF >
inline typename VT::ConstIterator begin( const Vector<VT,TF>& v )
{
   return (~v).begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given vector.
//
// \param v The given dense or sparse vector.
// \return Iterator to the first element of the given vector.
*/
template< typename VT, bool TF >
inline typename VT::ConstIterator cbegin( const Vector<VT,TF>& v )
{
   return (~v).begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given vector.
//
// \param v The given dense or sparse vector.
// \return Iterator just past the last element of the given vector.
*/
template< typename VT, bool TF >
inline typename VT::Iterator end( Vector<VT,TF>& v )
{
   return (~v).end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given vector.
//
// \param v The given dense or sparse vector.
// \return Iterator just past the last element of the given vector.
*/
template< typename VT, bool TF >
inline typename VT::ConstIterator end( const Vector<VT,TF>& v )
{
   return (~v).end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given vector.
//
// \param v The given dense or sparse vector.
// \return Iterator just past the last element of the given vector.
*/
template< typename VT, bool TF >
inline typename VT::ConstIterator cend( const Vector<VT,TF>& v )
{
   return (~v).end();
}
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
inline size_t size( const Vector<VT,TF>& v )
{
   return (~v).size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the vector.
// \ingroup vector
//
// \param v The given vector.
// \return The capacity of the vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
inline size_t capacity( const Vector<VT,TF>& v )
{
   return (~v).capacity();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the vector.
// \ingroup vector
//
// \param v The given vector.
// \return The number of non-zero elements in the vector.
//
// Note that the number of non-zero elements is always less than or equal to the current size
// of the vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
inline size_t nonZeros( const Vector<VT,TF>& v )
{
   return (~v).nonZeros();
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


//*************************************************************************************************
/*!\brief Returns whether the two given vectors represent the same observable state.
// \ingroup vector
//
// \param a The first vector to be tested for its state.
// \param b The second vector to be tested for its state.
// \return \a true in case the two vectors share a state, \a false otherwise.
//
// The isSame function provides an abstract interface for testing if the two given vectors
// represent the same observable state. This happens for instance in case \c a and \c b refer
// to the same vector or in case \c a and \c b are aliases for the same vector. In case both
// vectors represent the same observable state, the function returns \a true, other it returns
// \a false.

   \code
   typedef blaze::DynamicVector<int>          VectorType;
   typedef blaze::DenseSubvector<VectorType>  SubvectorType;

   VectorType vec1( 4UL );  // Setup of a 4-dimensional dynamic vector
   VectorType vec2( 4UL );  // Setup of a second 4-dimensional dynamic vector

   SubvectorType sub1 = subvector( vec1, 0UL, 4UL );  // Subvector of vec1 for the entire range
   SubvectorType sub2 = subvector( vec1, 1UL, 2UL );  // Subvector of vec1 for the range [1..3]
   SubvectorType sub3 = subvector( vec1, 1UL, 2UL );  // Second subvector of vec1 for the range [1..3]

   isSame( vec1, vec1 );  // returns true since both objects refer to the same vector
   isSame( vec1, vec2 );  // returns false since vec1 and vec2 are two different vectors
   isSame( vec1, sub1 );  // returns true since sub1 represents the same observable state as vec1
   isSame( vec1, sub3 );  // returns false since sub3 only covers part of the range of vec1
   isSame( sub2, sub3 );  // returns true since sub1 and sub2 refer to exactly the same range of vec1
   isSame( sub1, sub3 );  // returns false since sub1 and sub3 refer to different ranges of vec1
   \endcode
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline bool isSame( const Vector<VT1,TF1>& a, const Vector<VT2,TF2>& b )
{
   return ( IsSame<VT1,VT2>::value &&
            reinterpret_cast<const void*>( &a ) == reinterpret_cast<const void*>( &b ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
