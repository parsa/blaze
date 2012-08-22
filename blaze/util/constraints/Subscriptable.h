//=================================================================================================
/*!
//  \file blaze/util/constraints/Subscriptable.h
//  \brief Constraint on the data type
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
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
