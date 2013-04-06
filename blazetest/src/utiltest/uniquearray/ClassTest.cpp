//=================================================================================================
/*!
//  \file src/utiltest/uniquearray/ClassTest.cpp
//  \brief Source file for the UniqueArray class test
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <blaze/util/UniqueArray.h>
#include <blazetest/utiltest/uniquearray/ClassTest.h>
#include <blazetest/utiltest/Resource.h>


namespace blazetest {

namespace utiltest {

namespace uniquearray {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the UniqueArray class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testSingleResource();
   testRelease();
   testReset();
   testSelfReset();
   testSwap();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the general functionality of UniqueArray with a single resource.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs an initial general test of the UniqueArray functionality. It checks
// the number of Resource instances prior, during, and after passing a unique pointer a
// dynamically allocated Resource object. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSingleResource()
{
   // Initial check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Initial check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }

   // Acquiring a resource
   {
      blaze::UniqueArray<Resource> array( new Resource[3] );

      if( Resource::getCount() != 3U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 3\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Acquiring the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Final check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Final check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the release function of the UniqueArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the release function of the UniqueArray class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testRelease()
{
   // Initial check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Initial check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }

   {
      // Acquiring a resource
      blaze::UniqueArray<Resource> array( new Resource[4] );

      if( Resource::getCount() != 4U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 4\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Acquiring the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Releasing the resource
      Resource* resource = array.release();

      if( Resource::getCount() != 4U ) {
         std::ostringstream oss;
         oss << " Test: Releasing the resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 4\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array.get() != NULL ) {
         std::ostringstream oss;
         oss << " Test: Releasing the resource\n"
             << " Error: Releasing the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Manual destruction of the resource
      delete[] resource;
   }

   // Final check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Final check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reset function of the UniqueArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset function of the UniqueArray class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   // Initial check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Initial check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }

   {
      // Acquiring a resource
      blaze::UniqueArray<Resource> array( new Resource[5] );

      if( Resource::getCount() != 5U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 5\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Acquiring the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the resource
      array.reset();

      if( Resource::getCount() != 0U ) {
         std::ostringstream oss;
         oss << " Test: Resetting the resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 0\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array.get() != NULL ) {
         std::ostringstream oss;
         oss << " Test: Resetting the resource\n"
             << " Error: Resetting the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Final check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Final check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reset function of the UniqueArray class template with self assignment.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs an extended test of the release function of the UniqueArray class
// template involving self assignment. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSelfReset()
{
   // Initial check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Initial check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }

   {
      // Acquiring a resource
      blaze::UniqueArray<Resource> array( new Resource[6] );

      if( Resource::getCount() != 6U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 6\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Acquiring the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Self-resetting the unique array
      array.reset( array.get() );

      if( Resource::getCount() != 6U ) {
         std::ostringstream oss;
         oss << " Test: Self-resetting the unique array\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 6\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Self-resetting the resource\n"
             << " Error: Self-resetting the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Final check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Final check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the swap functionality of the UniqueArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the swap function of the UniqueArray class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   // Initial check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Initial check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }

   {
      // Acquiring two resources
      blaze::UniqueArray<Resource> array1( new Resource[3] );
      blaze::UniqueArray<Resource> array2( new Resource[5] );

      if( Resource::getCount() != 8U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring two resources\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 8\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array1.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring two resources\n"
             << " Error: Acquiring the resource for the first unique pointer failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array2.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring two resources\n"
             << " Error: Acquiring the resource for the second unique pointer failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Swapping the resources
      swap( array1, array2 );

      if( Resource::getCount() != 8U ) {
         std::ostringstream oss;
         oss << " Test: Swapping the resources\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 8\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array1.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Swapping the resources\n"
             << " Error: The first unique pointer was reset\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( array2.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Swapping the resources\n"
             << " Error: The second unique pointer was reset\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Final check of the resource counter
   if( Resource::getCount() != 0U ) {
      std::ostringstream oss;
      oss << " Test: Final check of the resource counter\n"
          << " Error: Invalid counter value\n"
          << " Details:\n"
          << "   Found counter    = " << Resource::getCount() << "\n"
          << "   Expected counter = 0\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************

} // namespace uniquearray

} // namespace utiltest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running UniqueArray class test..." << std::endl;

   try
   {
      RUN_UNIQUEARRAY_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UniqueArray class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
