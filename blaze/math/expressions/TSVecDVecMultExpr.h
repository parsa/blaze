//=================================================================================================
/*!
//  \file blaze/math/expressions/TSVecDVecMultExpr.h
//  \brief Header file for the sparse vector/dense vector inner product expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TSVECDVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TSVECDVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <boost/type_traits/remove_reference.hpp>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (inner product) of a sparse and a
//        dense vector (\f$ s=\vec{a}*\vec{b} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector for the inner product.
// \param rhs The right-hand side dense vector for the inner product.
// \return The scalar product.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the scalar product (inner product) of a sparse vector and a dense
// vector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::CompressedVector<double,rowVector> a;
   blaze::DynamicVector<double,columnVector> b;
   blaze::real res;
   // ... Resizing and initialization
   res = a * b;
   \endcode

// The operator returns a scalar value of the higher-order element type of the two involved
// vector element types \a T1::ElementType and \a T2::ElementType. Both vector types \a T1
// and \a T2 as well as the two element types \a T1::ElementType and \a T2::ElementType have
// to be supported by the MultTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1    // Type of the left-hand side sparse vector
        , typename T2 >  // Type of the right-hand side dense vector
inline const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type
   operator*( const SparseVector<T1,true>& lhs, const DenseVector<T2,false>& rhs )
{
   using boost::remove_reference;

   typedef typename T1::CompositeType            Lhs;            // Composite type of the left-hand side sparse vector expression
   typedef typename T2::CompositeType            Rhs;            // Composite type of the right-hand side dense vector expression
   typedef typename remove_reference<Lhs>::type  X1;             // Auxiliary type for the left-hand side composite type
   typedef typename remove_reference<Rhs>::type  X2;             // Auxiliary type for the right-hand side composite type
   typedef typename X1::ElementType              ET1;            // Element type of the left-hand side sparse vector expression
   typedef typename X2::ElementType              ET2;            // Element type of the right-hand side dense vector expression
   typedef typename MultTrait<ET1,ET2>::Type     MultType;       // Multiplication result type
   typedef typename X1::ConstIterator            ConstIterator;  // Iterator type of the left-hand sparse vector expression

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( T1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( T2 );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE( T1 );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( T2 );

   if( (~lhs).size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   if( (~lhs).nonZeros() == 0UL ) return MultType();

   Lhs left ( ~lhs );
   Rhs right( ~rhs );

   ConstIterator element( left.begin() );
   MultType sp( element->value() * right[ element->index() ] );
   ++element;

   for( ; element!=left.end(); ++element )
      sp += element->value() * right[ element->index() ];

   return sp;
}
//*************************************************************************************************

} // namespace blaze

#endif
