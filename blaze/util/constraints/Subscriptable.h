//=================================================================================================
/*!
//  \file blaze/util/constraints/Subscriptable.h
//  \brief Constraint on the data type
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

#ifndef _BLAZE_UTIL_CONSTRAINTS_SUBSCRIPTABLE_H_
#define _BLAZE_UTIL_CONSTRAINTS_SUBSCRIPTABLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/ConstraintTest.h>
#include <blaze/util/Suffix.h>


namespace blaze {

//=================================================================================================
//
//  MUST_BE_SUBSCRIBTABLE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the given type
// the \a constraint function can be called or not. If the given type is not subscriptable
// (the subscript operator can not be used on the type) a compilation error is created.
*/
template< typename T >
struct CONSTRAINT_MUST_BE_SUBSCRIBTABLE_FAILED
{
 private:
   //**********************************************************************************************
   static T createT();
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   enum { T_is_not_subscriptable = sizeof( createT()[0] ),
          value = 1 };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup constraints
//
// In case the given data type is not subscriptable, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_BE_SUBSCRIPTABLE(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_BE_SUBSCRIBTABLE_FAILED< T >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_BE_SUBSCRIPTABLE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_BE_SUBSCRIBTABLE_AS_DECAYABLE_POINTER CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the given type
// the \a constraint function can be called or not. If the given type is not a subscriptable
// pointer a compilation error is created.
*/
template< typename T >
struct CONSTRAINT_MUST_BE_SUBSCRIBTABLE_AS_DECAYABLE_POINTER_FAILED
{
 private:
   //**********************************************************************************************
   static T createT();
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   enum { T_is_not_subscriptable = sizeof( 0[createT()] ),
          value = 1 };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup constraints
//
// In case the given data type is not a subscriptable pointer, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_BE_SUBSCRIPTABLE_AS_DECAYABLE_POINTER(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_BE_SUBSCRIBTABLE_AS_DECAYABLE_POINTER_FAILED< T >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_BE_SUBSCRIBTABLE_AS_DECAYABLE_POINTER_TYPEDEF, __LINE__ )
//*************************************************************************************************

} // namespace blaze

#endif
