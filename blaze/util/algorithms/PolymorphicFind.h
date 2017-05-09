//=================================================================================================
/*!
//  \file blaze/util/algorithms/PolymorphicFind.h
//  \brief Headerfile for the generic polymorphicFind algorithm
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

#ifndef _BLAZE_UTIL_ALGORITHMS_POLYMORPHICFIND_H_
#define _BLAZE_UTIL_ALGORITHMS_POLYMORPHICFIND_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/DerivedFrom.h>


namespace blaze {

//=================================================================================================
//
//  POLYMORPHIC FIND ALGORITHM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Finds the next pointer to an object with dynamic type \a D.
// \ingroup algorithms
//
// \param first Iterator to the first pointer of the pointer range.
// \param last Iterator to the pointer one past the last pointer of the pointer range.
// \return The next pointer to an object with dynamic type \a D.
//
// This function traverses the range \f$ [first,last) \f$ of pointers to objects with static
// type \a S until it finds the next polymorphic pointer to an object of dynamic type \a D.
// Note that in case \a D is not a type derived from \a S, a compile time error is created!
*/
template< typename D    // Dynamic type of the objects
        , typename S >  // Static type of the objects
inline S *const * polymorphicFind( S *const * first, S *const * last )
{
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_DERIVED_FROM( D, S );

   while( first != last && !dynamic_cast<D*>( *first ) ) ++first;
   return first;
}
//*************************************************************************************************

} // namespace blaze

#endif
