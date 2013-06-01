//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/Complex.h
//  \brief Specialization of the Creator class template for complex values
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_COMPLEX_H_
#define _BLAZETEST_MATHTEST_CREATOR_COMPLEX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Boolean.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Void.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/creator/Default.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>
#include <blazetest/system/Types.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for complex data types.
//
// This specialization of the Creator class template is able to create random complex values.
*/
template< typename T >  // Element type of the complex type
class Creator< complex<T> >
{
 public:
   //**Type definitions****************************************************************************
   typedef complex<T>  Type;  //!< Type to be created by the Creator.
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
   const complex<T> operator()() const;
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
/*!\brief Returns a randomly created complex value.
//
// \return The randomly generated complex value.
*/
template< typename T >  // Element type of the complex type
inline const complex<T> Creator< complex<T> >::operator()() const
{
   return complex<T>( blaze::rand<T>( T(randmin), T(randmax) ) );
}
//*************************************************************************************************

} // namespace blazetest

#endif
