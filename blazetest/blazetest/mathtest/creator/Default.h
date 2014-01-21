//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/Default.h
//  \brief Header file for the Creator class template
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
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


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
// data type, a value in the range \f$ [0..1) \f$ is generated, in case \a T is a signed integral
// data type, the value will be in the range \f$ [-10..10] \f$, and in case \a T is an unsigned
// integral data type, a value in the range \f$ [0..10] \f$ is generated.
*/
template< typename T >  // Type to be created
inline T Creator<T>::operator()() const
{
   return blaze::rand<T>( T(randmin), T(randmax) );
}
//*************************************************************************************************

} // namespace blazetest

#endif
