//=================================================================================================
/*!
//  \file blaze/math/Band.h
//  \brief Header file for the complete Band implementation
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_BAND_H_
#define _BLAZE_MATH_BAND_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/Exception.h>
#include <blaze/math/smp/DenseVector.h>
#include <blaze/math/smp/SparseVector.h>
#include <blaze/math/views/Band.h>
#include <blaze/util/Random.h>


namespace blaze {

//=================================================================================================
//
//  RAND SPECIALIZATION FOR DENSE BANDS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for dense bands.
// \ingroup random
//
// This specialization of the Rand class randomizes dense bands.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
class Rand< DenseBand<MT,BIs...> >
{
 public:
   //**Randomize functions*************************************************************************
   /*!\name Randomize functions */
   //@{
   inline void randomize( DenseBand<MT,BIs...>& band ) const;

   template< typename Arg >
   inline void randomize( DenseBand<MT,BIs...>& band, const Arg& min, const Arg& max ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a dense band.
//
// \param band The band to be randomized.
// \return void
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
inline void Rand< DenseBand<MT,BIs...> >::randomize( DenseBand<MT,BIs...>& band ) const
{
   using blaze::randomize;

   for( size_t i=0UL; i<band.size(); ++i ) {
      randomize( band[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a dense band.
//
// \param band The band to be randomized.
// \param min The smallest possible value for a band element.
// \param max The largest possible value for a band element.
// \return void
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
template< typename Arg >      // Min/max argument type
inline void Rand< DenseBand<MT,BIs...> >::randomize( DenseBand<MT,BIs...>& band,
                                                     const Arg& min, const Arg& max ) const
{
   using blaze::randomize;

   for( size_t i=0UL; i<band.size(); ++i ) {
      randomize( band[i], min, max );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  RAND SPECIALIZATION FOR SPARSE BANDS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for sparse bands.
// \ingroup random
//
// This specialization of the Rand class randomizes sparse bands.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
class Rand< SparseBand<MT,BIs...> >
{
 public:
   //**Randomize functions*************************************************************************
   /*!\name Randomize functions */
   //@{
   inline void randomize( SparseBand<MT,BIs...>& band ) const;
   inline void randomize( SparseBand<MT,BIs...>& band, size_t nonzeros ) const;

   template< typename Arg >
   inline void randomize( SparseBand<MT,BIs...>& band, const Arg& min, const Arg& max ) const;

   template< typename Arg >
   inline void randomize( SparseBand<MT,BIs...>& band, size_t nonzeros, const Arg& min, const Arg& max ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a sparse band.
//
// \param band The band to be randomized.
// \return void
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
inline void Rand< SparseBand<MT,BIs...> >::randomize( SparseBand<MT,BIs...>& band ) const
{
   using ElementType = ElementType_< SparseBand<MT,BIs...> >;

   const size_t size( band.size() );

   if( size == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*size ) ) );

   band.reset();
   band.reserve( nonzeros );

   while( band.nonZeros() < nonzeros ) {
      band[ rand<size_t>( 0UL, size-1UL ) ] = rand<ElementType>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a sparse band.
//
// \param band The band to be randomized.
// \param nonzeros The number of non-zero elements of the random band.
// \return void
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
inline void Rand< SparseBand<MT,BIs...> >::randomize( SparseBand<MT,BIs...>& band, size_t nonzeros ) const
{
   using ElementType = ElementType_< SparseBand<MT,BIs...> >;

   const size_t size( band.size() );

   if( nonzeros > size ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
   }

   if( size == 0UL ) return;

   band.reset();
   band.reserve( nonzeros );

   while( band.nonZeros() < nonzeros ) {
      band[ rand<size_t>( 0UL, size-1UL ) ] = rand<ElementType>();
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a sparse band.
//
// \param band The band to be randomized.
// \param min The smallest possible value for a band element.
// \param max The largest possible value for a band element.
// \return void
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
template< typename Arg >      // Min/max argument type
inline void Rand< SparseBand<MT,BIs...> >::randomize( SparseBand<MT,BIs...>& band,
                                                      const Arg& min, const Arg& max ) const
{
   using ElementType = ElementType_< SparseBand<MT,BIs...> >;

   const size_t size( band.size() );

   if( size == 0UL ) return;

   const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*size ) ) );

   band.reset();
   band.reserve( nonzeros );

   while( band.nonZeros() < nonzeros ) {
      band[ rand<size_t>( 0UL, size-1UL ) ] = rand<ElementType>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a sparse band.
//
// \param band The band to be randomized.
// \param nonzeros The number of non-zero elements of the random band.
// \param min The smallest possible value for a band element.
// \param max The largest possible value for a band element.
// \return void
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
template< typename MT         // Type of the matrix
        , ptrdiff_t... BIs >  // Band indices
template< typename Arg >      // Min/max argument type
inline void Rand< SparseBand<MT,BIs...> >::randomize( SparseBand<MT,BIs...>& band, size_t nonzeros,
                                                      const Arg& min, const Arg& max ) const
{
   using ElementType = ElementType_< SparseBand<MT,BIs...> >;

   const size_t size( band.size() );

   if( nonzeros > size ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
   }

   if( size == 0UL ) return;

   band.reset();
   band.reserve( nonzeros );

   while( band.nonZeros() < nonzeros ) {
      band[ rand<size_t>( 0UL, size-1UL ) ] = rand<ElementType>( min, max );
   }
}
/*! \endcond */
//*************************************************************************************************


} // namespace blaze

#endif
