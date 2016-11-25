//=================================================================================================
/*!
//  \file blaze/math/Matrix.h
//  \brief Header file for all basic Matrix functionality
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

#ifndef _BLAZE_MATH_MATRIX_H_
#define _BLAZE_MATH_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <ostream>
#include <blaze/math/Aliases.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Matrix.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Matrix functions */
//@{
template< typename MT, bool SO >
inline auto trace( const Matrix<MT,SO>& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the trace of the given square matrix.
// \ingroup matrix
//
// \param m Reference to a constant matrix object.
// \return The trace of the matrix.
// \exception std::invalid_argument Invalid input matrix for trace computation.
//
// This function computes the trace of the given square matrix, i.e. sums the elements on its
// diagonal:

            \f[ trace(A) = a_{11} + a_{22} + ... + a_{nn} = \sum_{i=1}^{n} a_{ii} \f]

// In case the given matrix is not a square matrix a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline auto trace( const Matrix<MT,SO>& m )
{
   using ET = ElementType_<MT>;

   if( !isSquare( ~m ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix for trace computation" );
   }

   if( (~m).rows() == 0UL ) {
      return ET();
   }

   ET tmp( (~m)(0UL,0UL) );

   for( size_t i=1UL; i<(~m).rows(); ++i ) {
      tmp += (~m)(i,i);
   }

   return tmp;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Matrix operators */
//@{
template< typename MT, bool SO >
inline std::ostream& operator<<( std::ostream& os, const Matrix<MT,SO>& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for dense and sparse matrices.
// \ingroup matrix
//
// \param os Reference to the output stream.
// \param m Reference to a constant matrix object.
// \return Reference to the output stream.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline std::ostream& operator<<( std::ostream& os, const Matrix<MT,SO>& m )
{
   CompositeType_<MT> tmp( ~m );

   for( size_t i=0UL; i<tmp.rows(); ++i ) {
      os << "( ";
      for( size_t j=0UL; j<tmp.columns(); ++j ) {
         os << std::setw(12) << tmp(i,j) << " ";
      }
      os << ")\n";
   }

   return os;
}
//*************************************************************************************************

} // namespace blaze

#endif
