//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/CompressedHermitian.h
//  \brief Specialization of the Creator class template for HermitianMatrix<CompressedMatrix>
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_COMPRESSEDHERMITIAN_H_
#define _BLAZETEST_MATHTEST_CREATOR_COMPRESSEDHERMITIAN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/shims/Real.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/creator/Default.h>
#include <blazetest/system/Types.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for Hermitian compressed matrices.
//
// This specialization of the Creator class template is able to create random Hermitian compressed
// matrices.
*/
template< typename T     // Element type of the compressed matrix
        , bool SO        // Storage order of the compressed matrix
        , typename CP >  // Creation policy
class Creator< blaze::HermitianMatrix< blaze::CompressedMatrix<T,SO> >, CP >
{
 public:
   //**Type definitions****************************************************************************
   //! Type to be created by the Creator.
   typedef blaze::HermitianMatrix< blaze::CompressedMatrix<T,SO> >  Type;

   //! Creation policy for the built-in elements.
   typedef CP  Policy;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T,CP>& elementCreator = Creator<T,CP>() );
   explicit inline Creator( size_t n, size_t nonzeros,
                            const Creator<T,CP>& elementCreator = Creator<T,CP>() );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.
   const blaze::HermitianMatrix< blaze::CompressedMatrix<T,SO> > operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t n_;          //!< The number of rows and columns of the Hermitian compressed matrix.
   size_t nonzeros_;   //!< The number of non-zero elements in the Hermitian compressed matrix.
   Creator<T,CP> ec_;  //!< Creator for the elements of the Hermitian compressed matrix.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the creator specialization for HermitianMatrix<CompressedMatrix>.
//
// \param elementCreator The creator for the elements of the Hermitian compressed matrix.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename T     // Element type of the compressed matrix
        , bool SO        // Storage order of the compressed matrix
        , typename CP >  // Creation policy
inline Creator< blaze::HermitianMatrix< blaze::CompressedMatrix<T,SO> >, CP >::Creator( const Creator<T,CP>& elementCreator )
   : n_( 3UL )              // The number of rows and columns of the Hermitian compressed matrix
   , nonzeros_( 3UL )       // The total number of non-zero elements in the Hermitian compressed matrix
   , ec_( elementCreator )  // Creator for the elements of the Hermitian compressed matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for HermitianMatrix<CompressedMatrix>.
//
// \param n The number of rows and columns of the compressed matrix.
// \param nonzeros The number of non-zero elements in the compressed matrix.
// \param elementCreator The creator for the elements of the compressed matrix.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename T     // Element type of the compressed matrix
        , bool SO        // Storage order of the compressed matrix
        , typename CP >  // Creation policy
inline Creator< blaze::HermitianMatrix< blaze::CompressedMatrix<T,SO> >, CP >::Creator( size_t n, size_t nonzeros,
                                                                                        const Creator<T,CP>& elementCreator )
   : n_( n )                // The number of rows and columns of the Hermitian compressed matrix
   , nonzeros_( nonzeros )  // The total number of non-zero elements in the Hermitian compressed matrix
   , ec_( elementCreator )  // Creator for the elements of the Hermitian compressed matrix
{
   if( n_ * n_ < nonzeros_ )
      throw std::invalid_argument( "Invalid number of non-zero elements" );
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created Hermitian compressed matrix.
//
// \return The randomly generated Hermitian compressed matrix.
*/
template< typename T     // Element type of the compressed matrix
        , bool SO        // Storage order of the compressed matrix
        , typename CP >  // Creation policy
inline const blaze::HermitianMatrix< blaze::CompressedMatrix<T,SO> >
   Creator< blaze::HermitianMatrix< blaze::CompressedMatrix<T,SO> >, CP >::operator()() const
{
   using blaze::real;

   blaze::HermitianMatrix< blaze::CompressedMatrix<T,SO> > matrix( n_, nonzeros_ );

   while( matrix.nonZeros() < nonzeros_ )
   {
      const size_t row   ( blaze::rand<size_t>( 0UL, n_-1UL ) );
      const size_t column( blaze::rand<size_t>( 0UL, n_-1UL ) );

      if( row == column )
         matrix(row,column) = real( ec_() );
      else
         matrix(row,column) = ec_();
   }

   return matrix;
}
//*************************************************************************************************

} // namespace blazetest

#endif
