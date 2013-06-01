//=================================================================================================
/*!
//  \file blazetest/mathtest/RandomMaximum.h
//  \brief Header file for the RandomMaximum class template
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

#ifndef _BLAZETEST_MATHTEST_RANDOMMAXIMUM_H_
#define _BLAZETEST_MATHTEST_RANDOMMAXIMUM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <new>
#include <blazetest/mathtest/RandomLimits.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS TEMPLATE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Maximum value for random initializations.
//
// The RandomMaximum class is a wrapper class around the functionality of the RandomLimits class
// template. It represents the largest allowed value for a random initialization of a numeric
// data type. The RandomMaximum class can be implicitly converted to all built-in, numeric data
// types.\n
// In order to handle maximum random values conveniently, the global RandomMaximum instance
// blazetest::randmax is provided, which can be used wherever a numeric built-in data type is
// required:

   \code
   int    i = randmax;  // Assigns the largest allowed signed integer value
   float  f = randmax;  // Assigns the largest allowed single precision value
   double d = randmax;  // Assigns the largest allowed double precision value
   \endcode
*/
class RandomMaximum
{
 public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline RandomMaximum();
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
   RandomMaximum& operator=( const RandomMaximum& );
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
/*!\brief The default constructor of the RandomMaximum class.
*/
inline RandomMaximum::RandomMaximum()
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
// The conversion operator returns the largest allowed random value for the given data type \a T.
*/
template< typename T >
inline RandomMaximum::operator const T() const
{
   return RandomLimits<T>::max();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RANDOMMAXIMUM VALUE
//
//=================================================================================================

//*************************************************************************************************
const RandomMaximum randmax;
//*************************************************************************************************

} // namespace blazetest

#endif
