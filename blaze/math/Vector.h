//=================================================================================================
/*!
//  \file blaze/math/Vector.h
//  \brief Header file for all basic Vector functionality
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

#ifndef _BLAZE_MATH_VECTOR_H_
#define _BLAZE_MATH_VECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <ostream>
#include <blaze/math/expressions/Vector.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/TransposeFlag.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Vector operators */
//@{
template< typename T1, typename T2 >
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator,( const Vector<T1,false>& lhs, const Vector<T2,false>& rhs );

template< typename T1, typename T2 >
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator,( const Vector<T1,false>& lhs, const Vector<T2,true>& rhs );

template< typename T1, typename T2 >
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator,( const Vector<T1,true>& lhs, const Vector<T2,false>& rhs );

template< typename T1, typename T2 >
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator,( const Vector<T1,true>& lhs, const Vector<T2,true>& rhs );

template< typename VT, bool TF >
inline std::ostream& operator<<( std::ostream& os, const Vector<VT,TF>& dv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the inner product.
// \param rhs The right-hand side vector for the inner product.
// \return The scalar product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator,( const Vector<T1,false>& lhs, const Vector<T2,false>& rhs )
{
   return trans(~lhs) * (~rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the inner product.
// \param rhs The right-hand side vector for the inner product.
// \return The scalar product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator,( const Vector<T1,false>& lhs, const Vector<T2,true>& rhs )
{
   return trans(~lhs) * trans(~rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the inner product.
// \param rhs The right-hand side vector for the inner product.
// \return The scalar product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator,( const Vector<T1,true>& lhs, const Vector<T2,false>& rhs )
{
   return (~lhs) * (~rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the inner product.
// \param rhs The right-hand side vector for the inner product.
// \return The scalar product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator,( const Vector<T1,true>& lhs, const Vector<T2,true>& rhs )
{
   return (~lhs) * trans(~rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for dense and sparse vectors.
// \ingroup vector
//
// \param os Reference to the output stream.
// \param v Reference to a constant vector object.
// \return Reference to the output stream.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline std::ostream& operator<<( std::ostream& os, const Vector<VT,TF>& v )
{
   if( (~v).size() == 0UL ) {
      os << "( )\n";
   }
   else if( TF == rowVector ) {
      os << "(";
      for( size_t i=0UL; i<(~v).size(); ++i )
         os << " " << (~v)[i];
      os << " )\n";
   }
   else {
      for( size_t i=0UL; i<(~v).size(); ++i )
         os << "( " << std::setw( 11UL ) << (~v)[i] << " )\n";
   }

   return os;
}
//*************************************************************************************************

} // namespace blaze

#endif
