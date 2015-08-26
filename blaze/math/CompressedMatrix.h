//=================================================================================================
/*!
//  \file blaze/math/CompressedMatrix.h
//  \brief Header file for the complete CompressedMatrix implementation
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

#ifndef _BLAZE_MATH_COMPRESSEDMATRIX_H_
#define _BLAZE_MATH_COMPRESSEDMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <blaze/math/sparse/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/SparseMatrix.h>
#include <blaze/system/Precision.h>
#include <blaze/util/Exception.h>
#include <blaze/util/Random.h>


namespace blaze {

//=================================================================================================
//
//  RAND SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for CompressedMatrix.
// \ingroup random
//
// This specialization of the Rand class creates random instances of CompressedMatrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
class Rand< CompressedMatrix<Type,SO> >
{
 public:
   //**Generate functions**************************************************************************
   /*!\name Generate functions */
   //@{
   inline const CompressedMatrix<Type,SO> generate( size_t m, size_t n ) const;
   inline const CompressedMatrix<Type,SO> generate( size_t m, size_t n, size_t nonzeros ) const;

   template< typename Arg >
   inline const CompressedMatrix<Type,SO> generate( size_t m, size_t n, const Arg& min, const Arg& max ) const;

   template< typename Arg >
   inline const CompressedMatrix<Type,SO> generate( size_t m, size_t n, size_t nonzeros,
                                                    const Arg& min, const Arg& max ) const;
   //@}
   //**********************************************************************************************

   //**Randomize functions*************************************************************************
   /*!\name Randomize functions */
   //@{
   inline void randomize( CompressedMatrix<Type,SO>& matrix ) const;
   inline void randomize( CompressedMatrix<Type,SO>& matrix, size_t nonzeros ) const;

   template< typename Arg >
   inline void randomize( CompressedMatrix<Type,SO>& matrix, const Arg& min, const Arg& max ) const;

   template< typename Arg >
   inline void randomize( CompressedMatrix<Type,SO>& matrix, size_t nonzeros,
                          const Arg& min, const Arg& max ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random CompressedMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
// \return The generated random matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline const CompressedMatrix<Type,SO>
   Rand< CompressedMatrix<Type,SO> >::generate( size_t m, size_t n ) const
{
   CompressedMatrix<Type,SO> matrix( m, n );
   randomize( matrix );

   return matrix;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random CompressedMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
// \param nonzeros The number of non-zero elements of the random matrix.
// \return The generated random matrix.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline const CompressedMatrix<Type,SO>
   Rand< CompressedMatrix<Type,SO> >::generate( size_t m, size_t n, size_t nonzeros ) const
{
   if( nonzeros > m*n ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
   }

   CompressedMatrix<Type,SO> matrix( m, n );
   randomize( matrix, nonzeros );

   return matrix;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random CompressedMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
// \param min The smallest possible value for a matrix element.
// \return The generated random matrix.
// \param max The largest possible value for a matrix element.
*/
template< typename Type   // Data type of the matrix
        , bool SO >       // Storage order
template< typename Arg >  // Min/max argument type
inline const CompressedMatrix<Type,SO>
   Rand< CompressedMatrix<Type,SO> >::generate( size_t m, size_t n, const Arg& min, const Arg& max ) const
{
   CompressedMatrix<Type,SO> matrix( m, n );
   randomize( matrix, min, max );

   return matrix;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Generation of a random CompressedMatrix.
//
// \param m The number of rows of the random matrix.
// \param n The number of columns of the random matrix.
// \param nonzeros The number of non-zero elements of the random matrix.
// \param min The smallest possible value for a matrix element.
// \param max The largest possible value for a matrix element.
// \return The generated random matrix.
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type   // Data type of the matrix
        , bool SO >       // Storage order
template< typename Arg >  // Min/max argument type
inline const CompressedMatrix<Type,SO>
   Rand< CompressedMatrix<Type,SO> >::generate( size_t m, size_t n, size_t nonzeros,
                                                const Arg& min, const Arg& max ) const
{
   if( nonzeros > m*n ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
   }

   CompressedMatrix<Type,SO> matrix( m, n );
   randomize( matrix, nonzeros, min, max );

   return matrix;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a CompressedMatrix.
//
// \param matrix The matrix to be randomized.
// \return void
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void Rand< CompressedMatrix<Type,SO> >::randomize( CompressedMatrix<Type,SO>& matrix ) const
{
   const size_t m( matrix.rows()    );
   const size_t n( matrix.columns() );

   if( m == 0UL || n == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*m*n ) ) );

   matrix.reset();
   matrix.reserve( nonzeros );

   while( matrix.nonZeros() < nonzeros ) {
      matrix( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a CompressedMatrix.
//
// \param matrix The matrix to be randomized.
// \param nonzeros The number of non-zero elements of the random matrix.
// \return void
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void Rand< CompressedMatrix<Type,SO> >::randomize( CompressedMatrix<Type,SO>& matrix, size_t nonzeros ) const
{
   const size_t m( matrix.rows()    );
   const size_t n( matrix.columns() );

   if( nonzeros > m*n ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
   }

   if( m == 0UL || n == 0UL ) return;

   matrix.reset();
   matrix.reserve( nonzeros );

   while( matrix.nonZeros() < nonzeros ) {
      matrix( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<Type>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a CompressedMatrix.
//
// \param matrix The matrix to be randomized.
// \param min The smallest possible value for a matrix element.
// \param max The largest possible value for a matrix element.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO >       // Storage order
template< typename Arg >  // Min/max argument type
inline void Rand< CompressedMatrix<Type,SO> >::randomize( CompressedMatrix<Type,SO>& matrix,
                                                          const Arg& min, const Arg& max ) const
{
   const size_t m( matrix.rows()    );
   const size_t n( matrix.columns() );

   if( m == 0UL || n == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*m*n ) ) );

   matrix.reset();
   matrix.reserve( nonzeros );

   while( matrix.nonZeros() < nonzeros ) {
      matrix( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<Type>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a CompressedMatrix.
//
// \param matrix The matrix to be randomized.
// \param nonzeros The number of non-zero elements of the random matrix.
// \param min The smallest possible value for a matrix element.
// \param max The largest possible value for a matrix element.
// \return void
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename Type   // Data type of the matrix
        , bool SO >       // Storage order
template< typename Arg >  // Min/max argument type
inline void Rand< CompressedMatrix<Type,SO> >::randomize( CompressedMatrix<Type,SO>& matrix,
                                                          size_t nonzeros, const Arg& min, const Arg& max ) const
{
   const size_t m( matrix.rows()    );
   const size_t n( matrix.columns() );

   if( nonzeros > m*n ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
   }

   if( m == 0UL || n == 0UL ) return;

   matrix.reset();
   matrix.reserve( nonzeros );

   while( matrix.nonZeros() < nonzeros ) {
      matrix( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<Type>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief MxN single precision matrix.
// \ingroup compressed_matrix
*/
typedef CompressedMatrix<float,false>  CMatMxNf;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief MxN double precision matrix.
// \ingroup compressed_matrix
*/
typedef CompressedMatrix<double,false>  CMatMxNd;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief MxN matrix with system-specific precision.
// \ingroup compressed_matrix
*/
typedef CompressedMatrix<real_t,false>  CMatMxN;
//*************************************************************************************************

} // namespace blaze

#endif
