//=================================================================================================
/*!
//  \file blazetest/mathtest/RandomMinimum.h
//  \brief Header file for the RandomMinimum class template
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

#ifndef _BLAZETEST_MATHTEST_RANDOMMINIMUM_H_
#define _BLAZETEST_MATHTEST_RANDOMMINIMUM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <new>
#include <blazetest/mathtest/RandomLimits.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Minimum value for random initializations.
//
// The RandomMinimum class is a wrapper class around the functionality of the RandomLimits class
// template. It represents the smallest allowed value for a random initialization of a numeric
// data type. The RandomMinimum class can be implicitly converted to all built-in, numeric data
// types.\n
// In order to handle minimum random values conveniently, the global RandomMinimum instance
// blazetest::randmin is provided, which can be used wherever a numeric built-in data type is
// required:

   \code
   int    i = randmin;  // Assigns the smallest allowed signed integer value
   float  f = randmin;  // Assigns the smallest allowed single precision value
   double d = randmin;  // Assigns the smallest allowed double precision value
   \endcode
*/
class RandomMinimum
{
 public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline RandomMinimum();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   template< typename T >
   inline operator const T() const;
   //@}
   //**********************************************************************************************

 private:
   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   RandomMinimum& operator=( const RandomMinimum& );
   void* operator&() const;

   void* operator new  ( std::size_t ) /*throw( std::bad_alloc )*/;
   void* operator new[]( std::size_t ) /*throw( std::bad_alloc )*/;
   void* operator new  ( std::size_t, const std::nothrow_t& ) /*throw()*/;
   void* operator new[]( std::size_t, const std::nothrow_t& ) /*throw()*/;

   void operator delete  ( void* ) /*throw()*/;
   void operator delete[]( void* ) /*throw()*/;
   void operator delete  ( void*, const std::nothrow_t& ) /*throw()*/;
   void operator delete[]( void*, const std::nothrow_t& ) /*throw()*/;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor of the RandomMinimum class.
*/
inline RandomMinimum::RandomMinimum()
{}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion operator to the required data type.
//
// The conversion operator returns the smallest allowed random value for the given data type \a T.
*/
template< typename T >
inline RandomMinimum::operator const T() const
{
   return RandomLimits<T>::min();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RANDOMMINIMUM VALUE
//
//=================================================================================================

//*************************************************************************************************
const RandomMinimum randmin;
//*************************************************************************************************

} // namespace blazetest

#endif
