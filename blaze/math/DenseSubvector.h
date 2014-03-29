//=================================================================================================
/*!
//  \file blaze/math/DenseSubvector.h
//  \brief Header file for the complete DenseSubvector implementation
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

#ifndef _BLAZE_MATH_DENSESUBVECTOR_H_
#define _BLAZE_MATH_DENSESUBVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/smp/DenseVector.h>
#include <blaze/math/smp/SparseVector.h>
#include <blaze/math/views/DenseSubmatrix.h>
#include <blaze/math/views/DenseSubvector.h>
#include <blaze/math/views/SparseSubmatrix.h>
#include <blaze/math/views/SparseSubvector.h>
#include <blaze/math/views/Submatrix.h>
#include <blaze/math/views/Subvector.h>
#include <blaze/util/Random.h>


namespace blaze {

//=================================================================================================
//
//  RAND SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for DenseSubvector.
// \ingroup random
//
// This specialization of the Rand class randomizes instances of DenseSubvector.
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
class Rand< DenseSubvector<VT,AF,TF> >
{
 public:
   //**Randomize functions*************************************************************************
   /*!\name Randomize functions */
   //@{
   inline void randomize( DenseSubvector<VT,AF,TF>& subvector ) const;

   template< typename Arg >
   inline void randomize( DenseSubvector<VT,AF,TF>& subvector, const Arg& min, const Arg& max ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a DenseSubvector.
//
// \param subvector The subvector to be randomized.
// \return void
*/
template< typename VT  // Type of the dense vector
        , bool AF      // Alignment flag
        , bool TF >    // Transpose flag
inline void Rand< DenseSubvector<VT,AF,TF> >::randomize( DenseSubvector<VT,AF,TF>& subvector ) const
{
   using blaze::randomize;

   for( size_t i=0UL; i<subvector.size(); ++i ) {
      randomize( subvector[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a DenseSubvector.
//
// \param subvector The subvector to be randomized.
// \param min The smallest possible value for a subvector element.
// \param max The largest possible value for a subvector element.
// \return void
*/
template< typename VT     // Type of the dense vector
        , bool AF         // Alignment flag
        , bool TF >       // Transpose flag
template< typename Arg >  // Min/max argument type
inline void Rand< DenseSubvector<VT,AF,TF> >::randomize( DenseSubvector<VT,AF,TF>& subvector,
                                                         const Arg& min, const Arg& max ) const
{
   using blaze::randomize;

   for( size_t i=0UL; i<subvector.size(); ++i ) {
      randomize( subvector[i], min, max );
   }
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
