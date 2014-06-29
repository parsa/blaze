//=================================================================================================
/*!
//  \file blaze/math/DenseColumn.h
//  \brief Header file for the complete DenseColumn implementation
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

#ifndef _BLAZE_MATH_DENSECOLUMN_H_
#define _BLAZE_MATH_DENSECOLUMN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/smp/DenseVector.h>
#include <blaze/math/smp/SparseVector.h>
#include <blaze/math/views/Column.h>
#include <blaze/math/views/DenseColumn.h>
#include <blaze/math/views/DenseRow.h>
#include <blaze/math/views/Row.h>
#include <blaze/util/Random.h>


namespace blaze {

//=================================================================================================
//
//  RAND SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for DenseColumn.
// \ingroup random
//
// This specialization of the Rand class randomizes instances of DenseColumn.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
class Rand< DenseColumn<MT,SO> >
{
 public:
   //**Randomize functions*************************************************************************
   /*!\name Randomize functions */
   //@{
   inline void randomize( DenseColumn<MT,SO>& column ) const;

   template< typename Arg >
   inline void randomize( DenseColumn<MT,SO>& column, const Arg& min, const Arg& max ) const;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a DenseColumn.
//
// \param column The column to be randomized.
// \return void
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline void Rand< DenseColumn<MT,SO> >::randomize( DenseColumn<MT,SO>& column ) const
{
   using blaze::randomize;

   for( size_t i=0UL; i<column.size(); ++i ) {
      randomize( column[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Randomization of a DenseColumn.
//
// \param column The column to be randomized.
// \param min The smallest possible value for a column element.
// \param max The largest possible value for a column element.
// \return void
*/
template< typename MT     // Type of the dense matrix
        , bool SO >       // Storage order
template< typename Arg >  // Min/max argument type
inline void Rand< DenseColumn<MT,SO> >::randomize( DenseColumn<MT,SO>& column,
                                                   const Arg& min, const Arg& max ) const
{
   using blaze::randomize;

   for( size_t i=0UL; i<column.size(); ++i ) {
      randomize( column[i], min, max );
   }
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
