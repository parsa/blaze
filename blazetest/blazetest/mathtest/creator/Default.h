//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/Default.h
//  \brief Header file for the Creator class template
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_DEFAULT_H_
#define _BLAZETEST_MATHTEST_CREATOR_DEFAULT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Boolean.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Void.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default creator for random built-in data values.
//
// The Creator class creates random values of the given data type \a T. In case \a T is an
// integral data type, Creator returns values in the range \f$ [0..10] \f$, in case \a T is
// a floating point data type \f$ [0..1) \f$.
*/
template< typename T >  // Type to be created
class Creator
{
 public:
   //**Type definitions****************************************************************************
   typedef T  Type;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   // No explicitly declared constructor.
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.
   inline T operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE    ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST       ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE    ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOID        ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_BOOLEAN_TYPE( T );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created numeric value.
//
// \return The randomly generated numeric value.
//
// This operator returns a randomly generated numeric value. In case \a T is a floating point
// data type, a value in the range \f$ [0..1) \f$ is generated, in case \a T is an integral
// data type, a value in the range \f$ [0..10] \f$ is generated.
*/
template< typename T >  // Type to be created
inline T Creator<T>::operator()() const
{
   if( blaze::IsFloatingPoint<T>::value )
      return blaze::rand<T>();
   else
      return blaze::rand<T>( T(0), T(10) );
}
//*************************************************************************************************

} // namespace blazetest

#endif
