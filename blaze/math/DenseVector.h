//=================================================================================================
/*!
//  \file blaze/math/DenseVector.h
//  \brief Header file for all basic DenseVector functionality
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

#ifndef _BLAZE_MATH_DENSEVECTOR_H_
#define _BLAZE_MATH_DENSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_reference.hpp>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/DVecAbsExpr.h>
#include <blaze/math/expressions/DVecDVecAddExpr.h>
#include <blaze/math/expressions/DVecDVecCrossExpr.h>
#include <blaze/math/expressions/DVecDVecMultExpr.h>
#include <blaze/math/expressions/DVecDVecSubExpr.h>
#include <blaze/math/expressions/DVecEvalExpr.h>
#include <blaze/math/expressions/DVecScalarDivExpr.h>
#include <blaze/math/expressions/DVecScalarMultExpr.h>
#include <blaze/math/expressions/DVecSVecAddExpr.h>
#include <blaze/math/expressions/DVecSVecCrossExpr.h>
#include <blaze/math/expressions/DVecSVecSubExpr.h>
#include <blaze/math/expressions/DVecTransExpr.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/SVecDVecCrossExpr.h>
#include <blaze/math/expressions/SVecDVecSubExpr.h>
#include <blaze/math/expressions/TDVecDVecMultExpr.h>
#include <blaze/math/Functions.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/TransposeFlag.h>
#include <blaze/math/Vector.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseVector operators */
//@{
template< typename T1, bool TF1, typename T2, bool TF2 >
inline bool operator==( const DenseVector<T1,TF1>& lhs, const DenseVector<T2,TF2>& rhs );

template< typename T1, bool TF1, typename T2, bool TF2 >
inline bool operator==( const DenseVector<T1,TF1>& lhs, const SparseVector<T2,TF2>& rhs );

template< typename T1, bool TF1, typename T2, bool TF2 >
inline bool operator==( const SparseVector<T1,TF1>& lhs, const DenseVector<T2,TF2>& rhs );

template< typename T1, typename T2, bool TF >
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator==( const DenseVector<T1,TF>& vec, T2 scalar );

template< typename T1, typename T2, bool TF >
inline typename EnableIf< IsNumeric<T1>, bool >::Type
   operator==( T1 scalar, const DenseVector<T2,TF>& vec );

template< typename T1, bool TF1, typename T2, bool TF2 >
inline bool operator!=( const DenseVector<T1,TF1>& lhs, const DenseVector<T2,TF2>& rhs );

template< typename T1, bool TF1, typename T2, bool TF2 >
inline bool operator!=( const DenseVector<T1,TF1>& lhs, const SparseVector<T2,TF2>& rhs );

template< typename T1, bool TF1, typename T2, bool TF2 >
inline bool operator!=( const SparseVector<T1,TF1>& lhs, const DenseVector<T2,TF2>& rhs );

template< typename T1, typename T2, bool TF >
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator!=( const DenseVector<T1,TF>& vec, T2 scalar );

template< typename T1, typename T2, bool TF >
inline typename EnableIf< IsNumeric<T1>, bool >::Type
   operator!=( T1 scalar, const DenseVector<T2,TF>& vec );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two dense vectors.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the comparison.
// \param rhs The right-hand side dense vector for the comparison.
// \return \a true if the two vectors are equal, \a false if not.
*/
template< typename T1  // Type of the left-hand side dense vector
        , bool TF1     // Transpose flag of the left-hand side dense vector
        , typename T2  // Type of the right-hand side dense vector
        , bool TF2 >   // Transpose flag of the right-hand side dense vector
inline bool operator==( const DenseVector<T1,TF1>& lhs, const DenseVector<T2,TF2>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;

   // Early exit in case the vector sizes don't match
   if( (~lhs).size() != (~rhs).size() ) return false;

   // Evaluation of the two dense vector operands
   CT1 a( ~lhs );
   CT2 b( ~rhs );

   // In order to compare the two vectors, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t i=0; i<a.size(); ++i )
      if( !equal( a[i], b[i] ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a dense vector and a sparse vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the comparison.
// \param rhs The right-hand side sparse vector for the comparison.
// \return \a true if the two vectors are equal, \a false if not.
*/
template< typename T1  // Type of the left-hand side dense vector
        , bool TF1     // Transpose flag of the left-hand side dense vector
        , typename T2  // Type of the right-hand side sparse vector
        , bool TF2 >   // Transpose flag of the right-hand side sparse vector
inline bool operator==( const DenseVector<T1,TF1>& lhs, const SparseVector<T2,TF2>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;
   typedef typename boost::remove_reference<CT2>::type::ConstIterator  ConstIterator;

   // Early exit in case the vector sizes don't match
   if( (~lhs).size() != (~rhs).size() ) return false;

   // Evaluation of the dense vector and sparse vector operand
   CT1 a( ~lhs );
   CT2 b( ~rhs );

   // In order to compare the two vectors, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   size_t i( 0 );

   for( ConstIterator element=b.begin(); element!=b.end(); ++element, ++i ) {
      for( ; i<element->index(); ++i ) {
         if( !isDefault( a[i] ) ) return false;
      }
      if( !equal( element->value(), a[i] ) ) return false;
   }
   for( ; i<a.size(); ++i ) {
      if( !isDefault( a[i] ) ) return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a sparse vector and a dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side sparse vector for the comparison.
// \param rhs The right-hand side dense vector for the comparison.
// \return \a true if the two vectors are equal, \a false if not.
*/
template< typename T1  // Type of the left-hand side sparse vector
        , bool TF1     // Transpose flag of the left-hand side sparse vector
        , typename T2  // Type of the right-hand side dense vector
        , bool TF2 >   // Transpose flag of the right-hand side dense vector
inline bool operator==( const SparseVector<T1,TF1>& lhs, const DenseVector<T2,TF2>& rhs )
{
   return ( rhs == lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if all elements of the vector are equal to the scalar, \a false if not.
//
// If all values of the vector are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side scalar
        , bool TF >    // Transpose flag
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator==( const DenseVector<T1,TF>& vec, T2 scalar )
{
   typedef typename T1::CompositeType  CT1;

   // Evaluation of the dense vector operand
   CT1 a( ~vec );

   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   for( size_t i=0; i<a.size(); ++i )
      if( !equal( a[i], scalar ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a scalar value and a dense vector.
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the comparison.
// \param vec The right-hand side dense vector for the comparison.
// \return \a true if all elements of the vector are equal to the scalar, \a false if not.
//
// If all values of the vector are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsNumeric<T1>, bool >::Type
   operator==( T1 scalar, const DenseVector<T2,TF>& vec )
{
   return ( vec == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two dense vectors.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the comparison.
// \param rhs The right-hand side dense vector for the comparison.
// \return \a true if the two vectors are not equal, \a false if they are equal.
*/
template< typename T1  // Type of the left-hand side dense vector
        , bool TF1     // Transpose flag of the left-hand side dense vector
        , typename T2  // Type of the right-hand side dense vector
        , bool TF2 >   // Transpose flag of the right-hand side dense vector
inline bool operator!=( const DenseVector<T1,TF1>& lhs, const DenseVector<T2,TF2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a dense vector and a sparse vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the comparison.
// \param rhs The right-hand side sparse vector for the comparison.
// \return \a true if the two vectors are not equal, \a false if they are equal.
*/
template< typename T1  // Type of the left-hand side dense vector
        , bool TF1     // Transpose flag of the left-hand side dense vector
        , typename T2  // Type of the right-hand side sparse vector
        , bool TF2 >   // Transpose flag of the right-hand side sparse vector
inline bool operator!=( const DenseVector<T1,TF1>& lhs, const SparseVector<T2,TF2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a sparse vector and a dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side sparse vector for the comparison.
// \param rhs The right-hand side dense vector for the comparison.
// \return \a true if the two vectors are not equal, \a false if they are equal.
*/
template< typename T1  // Type of the left-hand side sparse vector
        , bool TF1     // Transpose flag of the left-hand side sparse vector
        , typename T2  // Type of the right-hand side dense vector
        , bool TF2 >   // Transpose flag of the right-hand side dense vector
inline bool operator!=( const SparseVector<T1,TF1>& lhs, const DenseVector<T2,TF2>& rhs )
{
   return !( rhs == lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if at least one element of the vector is different from the scalar, \a false if not.
//
// If one value of the vector is inequal to the scalar value, the inequality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side scalar
        , bool TF >    // Transpose flag
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator!=( const DenseVector<T1,TF>& vec, T2 scalar )
{
   return !( vec == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a scalar value and a dense vector.
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the comparison.
// \param vec The right-hand side dense vector for the comparison.
// \return \a true if at least one element of the vector is different from the scalar, \a false if not.
//
// If one value of the vector is inequal to the scalar value, the inequality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side vector
        , bool TF >    // Transpose flag
inline typename EnableIf< IsNumeric<T1>, bool >::Type
   operator!=( T1 scalar, const DenseVector<T2,TF>& vec )
{
   return !( vec == scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseVector functions */
//@{
template< typename VT, bool TF >
inline const typename VT::ElementType min( const DenseVector<VT,TF>& dv );

template< typename VT, bool TF >
inline const typename VT::ElementType max( const DenseVector<VT,TF>& dv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the smallest element of the dense vector.
// \ingroup dense_vector
//
// \return The smallest dense vector element.
//
// This function returns the smallest element of the given dense vector. This function can
// only be used for element types that support the smaller-than relationship. In case the
// vector currently has a size of 0, the returned value is the default value (e.g. 0 in case
// of fundamental data types).
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline const typename VT::ElementType min( const DenseVector<VT,TF>& dv )
{
   using blaze::min;

   typedef typename VT::ElementType  ET;

   if( (~dv).size() == 0 ) return ET();

   ET minimum( (~dv)[0] );
   for( size_t i=1; i<(~dv).size(); ++i )
      minimum = min( minimum, (~dv)[i] );
   return minimum;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the largest element of the dense vector.
// \ingroup dense_vector
//
// \return The largest dense vector element.
//
// This function returns the largest element of the given dense vector. This function can
// only be used for element types that support the smaller-than relationship. In case the
// vector currently has a size of 0, the returned value is the default value (e.g. 0 in case
// of fundamental data types).
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline const typename VT::ElementType max( const DenseVector<VT,TF>& dv )
{
   using blaze::max;

   typedef typename VT::ElementType  ET;

   if( (~dv).size() == 0 ) return ET();

   ET maximum( (~dv)[0] );
   for( size_t i=1; i<(~dv).size(); ++i )
      maximum = max( maximum, (~dv)[i] );
   return maximum;
}
//*************************************************************************************************

} // namespace blaze

#endif
