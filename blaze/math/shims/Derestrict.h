//=================================================================================================
/*!
//  \file blaze/math/shims/Derestrict.h
//  \brief Header file for the derestrict shim
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

#ifndef _BLAZE_MATH_SHIMS_DERESTRICT_H_
#define _BLAZE_MATH_SHIMS_DERESTRICT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/UnrestrictedType.h>
#include <blaze/util/constraints/Volatile.h>


namespace blaze {

//=================================================================================================
//
//  DERESTRICT SHIM
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Formal removal of all restrictions on the data access to the given argument.
// \ingroup math_shims
//
// \param a The value/object to be derestricted.
// \return Reference to the argument without access restrictions.
//
// The derestrict shim represents an abstract interface for the removal of all restrictions on
// the data access to the given argument. It returns a reference to the argument that does
// provide the same interface but does not have any restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename T >
inline typename UnrestrictedType<T>::Type& derestrict( T& a )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( T );
   return a;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Formal removal of all restrictions on the data access to the given argument.
// \ingroup math_shims
//
// \param a The value/object to be derestricted.
// \return Reference to the argument without access restrictions.
//
// The derestrict shim represents an abstract interface for the removal of all restrictions on
// the data access to the given argument. It returns a reference to the argument that does
// provide the same interface but does not have any restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename T >
inline const typename UnrestrictedType<T>::Type& derestrict( const T& a )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( T );
   return derestrict( const_cast<T&>( a ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
