//=================================================================================================
/*!
//  \file blaze/math/shims/Invert.h
//  \brief Header file for the invert shim
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

#ifndef _BLAZE_MATH_SHIMS_INVERT_H_
#define _BLAZE_MATH_SHIMS_INVERT_H_


namespace blaze {

//=================================================================================================
//
//  INVERT SHIMS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Inverting the given single precision value.
// \ingroup math_shims
//
// \param a The single precision value to be inverted.
// \return The inverse of the given value.
//
// The invert shim represents an abstract interface for inverting a value/object of any given
// data type. For single precision floating point values this results in \f$ \frac{1}{a} \f$.
*/
inline float inv( float a )
{
   return ( 1.0F / a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given double precision value.
// \ingroup math_shims
//
// \param a The double precision value to be inverted.
// \return The inverse of the given value.
//
// The invert shim represents an abstract interface for inverting a value/object of any given
// data type. For double precision floating point values this results in \f$ \frac{1}{a} \f$.
*/
inline double inv( double a )
{
   return ( 1.0 / a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given long double value.
// \ingroup math_shims
//
// \param a The long double value to be inverted.
// \return The inverse of the given value.
//
// The invert shim represents an abstract interface for inverting a value/object of any given
// data type. For long double floating point values this results in \f$ \frac{1}{a} \f$.
*/
inline long double inv( long double a )
{
   return ( 1.0L / a );
}
//*************************************************************************************************

} // namespace blaze

#endif
