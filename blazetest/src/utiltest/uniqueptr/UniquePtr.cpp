//=================================================================================================
/*!
//  \file src/utiltest/uniqueptr/UniquePtr.cpp
//  \brief Source file for the UniquePtr test
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
#include <blaze/util/UniquePtr.h>
#include <blazetest/utiltest/UniquePtr.h>
#include <blazetest/utiltest/Resource.h>


namespace blazetest {

namespace utiltest {

namespace uniqueptr {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the UniquePtr class.
//
// \exception std::runtime_error Operation error detected.
*/
UniquePtr::UniquePtr()
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
/*!\brief Test of the general functionality of UniquePtr with a single resource.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs an initial general test of the UniquePtr functionality. It checks
// the number of Resource instances prior, during, and after passing a unique pointer a
// dynamically allocated Resource object. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void UniquePtr::testSingleResource()
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
      blaze::UniquePtr<Resource> ptr( new Resource() );

      if( Resource::getCount() != 1U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 1\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr.get() == NULL ) {
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
/*!\brief Test of the release function of the UniquePtr class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the release function of the UniquePtr class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UniquePtr::testRelease()
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
      blaze::UniquePtr<Resource> ptr( new Resource() );

      if( Resource::getCount() != 1U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 1\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Acquiring the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Releasing the resource
      Resource* resource = ptr.release();

      if( Resource::getCount() != 1U ) {
         std::ostringstream oss;
         oss << " Test: Releasing the resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 1\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr.get() != NULL ) {
         std::ostringstream oss;
         oss << " Test: Releasing the resource\n"
             << " Error: Releasing the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Manual destruction of the resource
      delete resource;
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
/*!\brief Test of the reset function of the UniquePtr class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset function of the UniquePtr class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UniquePtr::testReset()
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
      blaze::UniquePtr<Resource> ptr( new Resource() );

      if( Resource::getCount() != 1U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 1\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Acquiring the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the resource
      ptr.reset();

      if( Resource::getCount() != 0U ) {
         std::ostringstream oss;
         oss << " Test: Resetting the resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 0\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr.get() != NULL ) {
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
/*!\brief Test of the reset function of the UniquePtr class template with self assignment.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs an extended test of the release function of the UniquePtr class
// template involving self assignment. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void UniquePtr::testSelfReset()
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
      blaze::UniquePtr<Resource> ptr( new Resource() );

      if( Resource::getCount() != 1U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 1\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring a resource\n"
             << " Error: Acquiring the resource failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Self-resetting the unique ptr
      ptr.reset( ptr.get() );

      if( Resource::getCount() != 1U ) {
         std::ostringstream oss;
         oss << " Test: Self-resetting the unique ptr\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 1\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr.get() == NULL ) {
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
/*!\brief Test of the swap functionality of the UniquePtr class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the swap function of the UniquePtr class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UniquePtr::testSwap()
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
      blaze::UniquePtr<Resource> ptr1( new Resource() );
      blaze::UniquePtr<Resource> ptr2( new Resource() );

      if( Resource::getCount() != 2U ) {
         std::ostringstream oss;
         oss << " Test: Acquiring two resources\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 2\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr1.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring two resources\n"
             << " Error: Acquiring the resource for the first unique pointer failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr2.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Acquiring two resources\n"
             << " Error: Acquiring the resource for the second unique pointer failed\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Swapping the resources
      swap( ptr1, ptr2 );

      if( Resource::getCount() != 2U ) {
         std::ostringstream oss;
         oss << " Test: Swapping the resources\n"
             << " Error: Invalid counter value\n"
             << " Details:\n"
             << "   Found counter    = " << Resource::getCount() << "\n"
             << "   Expected counter = 2\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr1.get() == NULL ) {
         std::ostringstream oss;
         oss << " Test: Swapping the resources\n"
             << " Error: The first unique pointer was reset\n"
             << " Details:\n"
             << "   Instance counter = " << Resource::getCount() << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( ptr2.get() == NULL ) {
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

} // namespace uniqueptr

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
   std::cout << "   Running UniquePtr test..." << std::endl;

   try
   {
      RUN_UNIQUEPTR_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UniquePtr test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
