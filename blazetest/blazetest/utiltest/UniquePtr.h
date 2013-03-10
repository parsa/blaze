//=================================================================================================
/*!
//  \file blazetest/utiltest/UniquePtr.h
//  \brief Header file for the UniquePtr test
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

#ifndef _BLAZETEST_UTILTEST_UNIQUEPTR_H_
#define _BLAZETEST_UTILTEST_UNIQUEPTR_H_


namespace blazetest {

namespace utiltest {

namespace uniqueptr {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for the test of the UniquePtr class template.
//
// The UniquePtr class template represents the collection of tests for the UniquePtr class
// template.
*/
class UniquePtr
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit UniquePtr();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   void testSingleResource();
   void testRelease();
   void testReset();
   void testSelfReset();
   void testSwap();
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the functionality of the UniquePtr class template.
//
// \return void
*/
inline void runTest()
{
   UniquePtr();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the UniquePtr test.
*/
#define RUN_UNIQUEPTR_TEST \
   blazetest::utiltest::uniqueptr::runTest();
/*! \endcond */
//*************************************************************************************************

} // namespace uniqueptr

} // namespace utiltest

} // namespace blazetest

#endif
