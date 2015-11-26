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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/shims/Square.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Complex.h>


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
// The \a inv shim represents an abstract interface for inverting a value/object of any given
// data type. For single precision floating point values this results in \f$ \frac{1}{a} \f$.
*/
BLAZE_ALWAYS_INLINE float inv( float a )
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
// The \a inv shim represents an abstract interface for inverting a value/object of any given
// data type. For double precision floating point values this results in \f$ \frac{1}{a} \f$.
*/
BLAZE_ALWAYS_INLINE double inv( double a )
{
   return ( 1.0 / a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given extended precision value.
// \ingroup math_shims
//
// \param a The extended precision value to be inverted.
// \return The inverse of the given value.
//
// The \a inv shim represents an abstract interface for inverting a value/object of any given
// data type. For extended precision floating point values this results in \f$ \frac{1}{a} \f$.
*/
BLAZE_ALWAYS_INLINE long double inv( long double a )
{
   return ( 1.0L / a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given single precision complex number.
// \ingroup math_shims
//
// \param a The single precision complex number to be inverted.
// \return The inverse of the given value.
//
// The \a inv shim represents an abstract interface for inverting a value/object of any given
// data type. For a single precision floating point complex number \f$ z = x + yi \f$ this
// results in \f$ \frac{\overline{z}}{x^2+y^2} \f$.
*/
BLAZE_ALWAYS_INLINE complex<float> inv( const complex<float>& a )
{
   const float tmp( 1.0F / ( sq( real(a) ) + sq( imag(a) ) ) );
   return complex<float>( tmp*real(a), -tmp*imag(a) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given double precision complex number.
// \ingroup math_shims
//
// \param a The double precision complex number to be inverted.
// \return The inverse of the given value.
//
// The \a inv shim represents an abstract interface for inverting a value/object of any given
// data type. For a double precision floating point complex number \f$ z = x + yi \f$ this
// results in \f$ \frac{\overline{z}}{x^2+y^2} \f$.
*/
BLAZE_ALWAYS_INLINE complex<double> inv( const complex<double>& a )
{
   const double tmp( 1.0 / ( sq( real(a) ) + sq( imag(a) ) ) );
   return complex<double>( tmp*real(a), -tmp*imag(a) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given extended precision complex number.
// \ingroup math_shims
//
// \param a The extended precision complex number to be inverted.
// \return The inverse of the given value.
//
// The \a inv shim represents an abstract interface for inverting a value/object of any given
// data type. For an extended precision floating point complex number \f$ z = x + yi \f$ this
// results in \f$ \frac{\overline{z}}{x^2+y^2} \f$.
*/
BLAZE_ALWAYS_INLINE complex<long double> inv( const complex<long double>& a )
{
   const long double tmp( 1.0L / ( sq( real(a) ) + sq( imag(a) ) ) );
   return complex<long double>( tmp*real(a), -tmp*imag(a) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given single precision value.
// \ingroup math_shims
//
// \param a The single precision value to be inverted.
// \return The inverse of the given value.
//
// The \a invert shim represents an abstract interface for inverting a value/object of any
// given data type in-place. For single precision floating point values this results in
// \f$ \frac{1}{a} \f$.
*/
BLAZE_ALWAYS_INLINE void invert( float& a )
{
   a = inv( a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given double precision value.
// \ingroup math_shims
//
// \param a The double precision value to be inverted.
// \return The inverse of the given value.
//
// The \a invert shim represents an abstract interface for inverting a value/object of any
// given data type in-place. For double precision floating point values this results in
// \f$ \frac{1}{a} \f$.
*/
BLAZE_ALWAYS_INLINE void invert( double& a )
{
   a = inv( a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given extended precision value.
// \ingroup math_shims
//
// \param a The extended precision value to be inverted.
// \return The inverse of the given value.
//
// The \a invert shim represents an abstract interface for inverting a value/object of any
// given data type in-place. For extended precision floating point values this results in
// \f$ \frac{1}{a} \f$.
*/
BLAZE_ALWAYS_INLINE void invert( long double& a )
{
   a = inv( a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given single precision complex number.
// \ingroup math_shims
//
// \param a The single precision complex number to be inverted.
// \return The inverse of the given value.
//
// The \a invert shim represents an abstract interface for inverting a value/object of any given
// data type in-place. For a single precision floating point complex number \f$ z = x + yi \f$
// this results in \f$ \frac{\overline{z}}{x^2+y^2} \f$.
*/
BLAZE_ALWAYS_INLINE void invert( complex<float>& a )
{
   a = inv( a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given double precision complex number.
// \ingroup math_shims
//
// \param a The double precision complex number to be inverted.
// \return The inverse of the given value.
//
// The \a invert shim represents an abstract interface for inverting a value/object of any given
// data type in-place. For a double precision floating point complex number \f$ z = x + yi \f$
// this results in \f$ \frac{\overline{z}}{x^2+y^2} \f$.
*/
BLAZE_ALWAYS_INLINE void invert( complex<double>& a )
{
   a = inv( a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given extended precision complex number.
// \ingroup math_shims
//
// \param a The extended precision complex number to be inverted.
// \return The inverse of the given value.
//
// The \a invert shim represents an abstract interface for inverting a value/object of any given
// data type in-place. For an extended precision floating point complex number \f$ z = x + yi \f$
// this results in \f$ \frac{\overline{z}}{x^2+y^2} \f$.
*/
BLAZE_ALWAYS_INLINE void invert( complex<long double>& a )
{
   a = inv( a );
}
//*************************************************************************************************

} // namespace blaze

#endif
