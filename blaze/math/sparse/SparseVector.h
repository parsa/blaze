//=================================================================================================
/*!
//  \file blaze/math/sparse/SparseVector.h
//  \brief Header file for utility functions for sparse vectors
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_SPARSE_SPARSEVECTOR_H_
#define _BLAZE_MATH_SPARSE_SPARSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/shims/Pow2.h>
#include <blaze/math/shims/Sqrt.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingNumeric.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DecltypeAuto.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseVector operators */
//@{
template< typename VT, bool TF, typename ST >
inline auto operator*=( SparseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsNumeric_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
inline auto operator*=( SparseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsNumeric_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
inline auto operator/=( SparseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsNumeric_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
inline auto operator/=( SparseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsNumeric_v<ST>, VT& >;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a sparse vector and
//        a scalar value (\f$ \vec{a}*=s \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side sparse vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side sparse vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( SparseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsNumeric_v<ST>, VT& >
{
   if( IsRestricted_v<VT> ) {
      if( !tryMult( ~vec, 0UL, (~vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted vector" );
      }
   }

   if( !IsResizable_v< ElementType_t<VT> > && isZero( scalar ) )
   {
      reset( ~vec );
   }
   else
   {
      BLAZE_DECLTYPE_AUTO( left, derestrict( ~vec ) );

      const auto last( left.end() );
      for( auto element=left.begin(); element!=last; ++element ) {
         element->value() *= scalar;
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact( ~vec ), "Invariant violation detected" );

   return ~vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a temporary sparse vector
//        and a scalar (\f$ v*=s \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side temporary sparse vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side sparse vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( SparseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsNumeric_v<ST>, VT& >
{
   return operator*=( ~vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a sparse vector by a scalar value
//        (\f$ \vec{a}/=s \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side sparse vector for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side sparse vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( SparseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsNumeric_v<ST>, VT& >
{
   BLAZE_USER_ASSERT( !isZero( scalar ), "Division by zero detected" );

   if( IsRestricted_v<VT> ) {
      if( !tryDiv( ~vec, 0UL, (~vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted vector" );
      }
   }

   using ScalarType = If_t< IsFloatingPoint_v< UnderlyingBuiltin_t<VT> > ||
                            IsFloatingPoint_v< UnderlyingBuiltin_t<ST> >
                          , If_t< IsComplex_v< UnderlyingNumeric_t<VT> > && IsBuiltin_v<ST>
                                , DivTrait_t< UnderlyingBuiltin_t<VT>, ST >
                                , DivTrait_t< UnderlyingNumeric_t<VT>, ST > >
                          , ST >;

   BLAZE_DECLTYPE_AUTO( left, derestrict( ~vec ) );

   if( IsInvertible_v<ScalarType> ) {
      const ScalarType tmp( ScalarType(1)/static_cast<ScalarType>( scalar ) );
      for( auto element=left.begin(); element!=left.end(); ++element ) {
         element->value() *= tmp;
      }
   }
   else {
      for( auto element=left.begin(); element!=left.end(); ++element ) {
         element->value() /= scalar;
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact( ~vec ), "Invariant violation detected" );

   return ~vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a temporary sparse vector by a scalar
//        value (\f$ \vec{a}/=s \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side temporary sparse vector for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side sparse vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( SparseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsNumeric_v<ST>, VT& >
{
   return operator/=( ~vec, scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseVector functions */
//@{
template< typename VT, bool TF >
bool isnan( const SparseVector<VT,TF>& sv );

template< typename VT, bool TF >
bool isUniform( const SparseVector<VT,TF>& dv );

template< typename VT, bool TF >
const ElementType_t<VT> sqrLength( const SparseVector<VT,TF>& sv );

template< typename VT, bool TF >
inline auto length( const SparseVector<VT,TF>& sv ) -> decltype( sqrt( sqrLength( ~sv ) ) );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse vector for not-a-number elements.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be checked for not-a-number elements.
// \return \a true if at least one element of the vector is not-a-number, \a false otherwise.
//
// This function checks the N-dimensional sparse vector for not-a-number (NaN) elements. If
// at least one element of the vector is not-a-number, the function returns \a true, otherwise
// it returns \a false.

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isnan( a ) ) { ... }
   \endcode

// Note that this function only works for vectors with floating point elements. The attempt to
// use it for a vector with a non-floating point element type results in a compile time error.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline bool isnan( const SparseVector<VT,TF>& sv )
{
   using CT = CompositeType_t<VT>;
   using ConstIterator = ConstIterator_t< RemoveReference_t<CT> >;

   CT a( ~sv );  // Evaluation of the sparse vector operand

   const ConstIterator end( a.end() );
   for( ConstIterator element=a.begin(); element!=end; ++element ) {
      if( isnan( element->value() ) ) return true;
   }
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse vector is a uniform vector.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be checked.
// \return \a true if the vector is a uniform vector, \a false if not.
//
// This function checks if the given sparse vector is a uniform vector. The vector is considered
// to be uniform if all its elements are identical. The following code example demonstrates the
// use of the function:

   \code
   blaze::CompressedVector<int,blaze::columnVector> a, b;
   // ... Initialization
   if( isUniform( a ) ) { ... }
   \endcode

// It is also possible to check if a vector expression results in a uniform vector:

   \code
   if( isUniform( a + b ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary vector.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
bool isUniform( const SparseVector<VT,TF>& sv )
{
   using CT = CompositeType_t<VT>;
   using ConstReference = ConstReference_t< RemoveReference_t<CT> >;
   using ConstIterator = ConstIterator_t< RemoveReference_t<CT> >;

   if( IsUniform_v<VT> || (~sv).size() < 2UL )
      return true;

   CT a( ~sv );  // Evaluation of the sparse vector operand

   if( a.nonZeros() != a.size() )
   {
      for( ConstIterator element=a.begin(); element!=a.end(); ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }
   else
   {
      ConstReference cmp( a[0] );
      ConstIterator element( a.begin() );

      ++element;

      for( ; element!=a.end(); ++element ) {
         if( element->value() != cmp )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the square length (magnitude) of the sparse vector \f$|\vec{a}|^2\f$.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \return The square length (magnitude) of the vector.
//
// This function calculates the actual square length (magnitude) of the sparse vector.
//
// \note This operation is only defined for numeric data types. In case the element type is
// not a numeric data type (i.e. a user defined data type or boolean) the attempt to use the
// sqrLength() function results in a compile time error!
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
const ElementType_t<VT> sqrLength( const SparseVector<VT,TF>& sv )
{
   using ElementType   = ElementType_t<VT>;
   using ConstIterator = ConstIterator_t<VT>;

   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( ElementType );

   ElementType sum( 0 );
   for( ConstIterator element=(~sv).begin(); element!=(~sv).end(); ++element )
      sum += pow2( element->value() );
   return sum;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the length (magnitude) of the sparse vector \f$|\vec{a}|\f$.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \return The length (magnitude) of the sparse vector.
//
// This function calculates the actual length (magnitude) of the sparse vector. The return type
// of the length() function depends on the actual element type of the vector instance:
//
// <table border="0" cellspacing="0" cellpadding="1">
//    <tr>
//       <td width="250px"> \b Type </td>
//       <td width="100px"> \b LengthType </td>
//    </tr>
//    <tr>
//       <td>float</td>
//       <td>float</td>
//    </tr>
//    <tr>
//       <td>integral data types and double</td>
//       <td>double</td>
//    </tr>
//    <tr>
//       <td>long double</td>
//       <td>long double</td>
//    </tr>
//    <tr>
//       <td>complex<T></td>
//       <td>complex<T></td>
//    </tr>
// </table>
//
// \note This operation is only defined for numeric data types. In case the element type is
// not a numeric data type (i.e. a user defined data type or boolean) the attempt to use the
// length() function results in a compile time error!
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline auto length( const SparseVector<VT,TF>& sv ) -> decltype( sqrt( sqrLength( ~sv ) ) )
{
   return sqrt( sqrLength( ~sv ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
