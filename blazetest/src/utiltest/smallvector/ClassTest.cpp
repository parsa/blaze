//=================================================================================================
/*!
//  \file src/utiltest/smallvector/ClassTest.cpp
//  \brief Source file for the SmallVector class test
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <list>
#include <blaze/util/Random.h>
#include <blazetest/utiltest/smallvector/ClassTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace utiltest {

namespace smallvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SmallVector class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testConstructors();
   testAssignment();
   testSubscript();
   testAt();
   testIterator();
   testClear();
   testResize();
   testReserve();
   testShrinkToFit();
   testPushBack();
   testInsert();
   testErase();
   testSwap();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the SmallVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the SmallVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "SmallVector default constructor";

      blaze::SmallVector<int,5UL> vec;

      checkSize( vec, 0UL );
   }


   //=====================================================================================
   // Size constructor
   //=====================================================================================

   {
      test_ = "SmallVector size constructor (size 0)";

      blaze::SmallVector<int,5UL> vec( 0UL );

      checkSize( vec, 0UL );
   }

   {
      test_ = "SmallVector size constructor (size 4)";

      blaze::SmallVector<int,5UL> vec( 4UL );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
   }

   {
      test_ = "SmallVector size constructor (size 5)";

      blaze::SmallVector<int,5UL> vec( 5UL );

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
   }

   {
      test_ = "SmallVector size constructor (size 6)";

      blaze::SmallVector<int,5UL> vec( 6UL );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );
   }


   //=====================================================================================
   // Homogeneous initialization
   //=====================================================================================

   {
      test_ = "SmallVector homogeneous initialization constructor (size 0)";

      blaze::SmallVector<int,5UL> vec( 0UL, 2 );

      checkSize( vec, 0UL );
   }

   {
      test_ = "SmallVector homogeneous initialization constructor (size 4)";

      blaze::SmallVector<int,5UL> vec( 4UL, 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector homogeneous initialization constructor (size 5)";

      blaze::SmallVector<int,5UL> vec( 5UL, 2 );

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 || vec[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector homogeneous initialization constructor (size 6)";

      blaze::SmallVector<int,5UL> vec( 6UL, 2 );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 || vec[4] != 2 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Range initialization
   //=====================================================================================

   {
      test_ = "SmallVector range constructor (size 4)";

      std::list<int> list{ 1, 2, 3, 4 };
      blaze::SmallVector<int,5UL> vec( list.begin(), list.end() );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector range constructor (size 5)";

      std::list<int> list{ 1, 2, 3, 4, 5 };
      blaze::SmallVector<int,5UL> vec( list.begin(), list.end() );

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector range constructor (size 6)";

      std::list<int> list{ 1, 2, 3, 4, 5, 6 };
      blaze::SmallVector<int,6UL> vec( list.begin(), list.end() );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // List initialization
   //=====================================================================================

   {
      test_ = "SmallVector initializer list constructor (size 4)";

      blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector initializer list constructor (size 5)";

      blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5 };

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector initializer list constructor (size 6)";

      blaze::SmallVector<int,6UL> vec{ 1, 2, 3, 4, 5, 6 };

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "SmallVector copy constructor (size 0)";

      blaze::SmallVector<int,5UL> vec1( 0UL );
      blaze::SmallVector<int,5UL> vec2( vec1 );

      checkSize    ( vec2, 0UL );
      checkCapacity( vec2, 0UL );
   }

   {
      test_ = "SmallVector copy constructor (size 4)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4 };
      blaze::SmallVector<int,5UL> vec2( vec1 );

      checkSize    ( vec2, 4UL );
      checkCapacity( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector copy constructor (size 5)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5 };
      blaze::SmallVector<int,5UL> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector copy constructor (size 6)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5, 6 };
      blaze::SmallVector<int,5UL> vec2( vec1 );

      checkSize    ( vec2, 6UL );
      checkCapacity( vec2, 6UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 || vec2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Move constructor
   //=====================================================================================

   {
      test_ = "SmallVector move constructor (size 0)";

      blaze::SmallVector<int,5UL> vec1( 0UL );
      blaze::SmallVector<int,5UL> vec2( std::move( vec1 ) );

      checkSize    ( vec2, 0UL );
      checkCapacity( vec2, 0UL );
   }

   {
      test_ = "SmallVector move constructor (size 4)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4 };
      blaze::SmallVector<int,5UL> vec2( std::move( vec1 ) );

      checkSize    ( vec2, 4UL );
      checkCapacity( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector move constructor (size 5)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5 };
      blaze::SmallVector<int,5UL> vec2( std::move( vec1 ) );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector move constructor (size 6)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5, 6 };
      blaze::SmallVector<int,5UL> vec2( std::move( vec1 ) );

      checkSize    ( vec2, 6UL );
      checkCapacity( vec2, 6UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 || vec2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SmallVector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the SmallVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // List assignment
   //=====================================================================================

   {
      test_ = "SmallVector initializer list assignment (size 3 to 4)";

      blaze::SmallVector<int,5UL> vec{ 11, 12, 13 };
      vec = { 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector initializer list assignment (size 8 to 4)";

      blaze::SmallVector<int,5UL> vec{ 11, 12, 13, 14, 15, 16, 17, 18 };
      vec = { 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector initializer list assignment (size 3 to 5)";

      blaze::SmallVector<int,5UL> vec{ 11, 12, 13 };
      vec = { 1, 2, 3, 4, 5 };

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector initializer list assignment (size 8 to 5)";

      blaze::SmallVector<int,5UL> vec{ 11, 12, 13, 14, 15, 16, 17, 18 };
      vec = { 1, 2, 3, 4, 5 };

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector initializer list assignment (size 3 to 6)";

      blaze::SmallVector<int,5UL> vec{ 11, 12, 13 };
      vec = { 1, 2, 3, 4, 5, 6 };

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector initializer list assignment (size 8 to 6)";

      blaze::SmallVector<int,5UL> vec{ 11, 12, 13, 14, 15, 16, 17, 18 };
      vec = { 1, 2, 3, 4, 5, 6 };

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "SmallVector copy assignment (size 4)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4 };
      blaze::SmallVector<int,5UL> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 4UL );
      checkCapacity( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector copy assignment (size 5)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5 };
      blaze::SmallVector<int,5UL> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector copy assignment (size 6)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5, 6 };
      blaze::SmallVector<int,5UL> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 6UL );
      checkCapacity( vec2, 6UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 || vec2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector copy assignment stress test";

      blaze::SmallVector<int,5UL> vec1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 10UL ) );
         blaze::SmallVector<int,5UL> vec2( size );
         for( int& element : vec2 ) {
            element = blaze::rand<int>( min, max );
         }

         vec1 = vec2;

         if( vec1 != vec2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << vec1 << "\n"
                << "   Expected result:\n" << vec2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Move assignment
   //=====================================================================================

   {
      test_ = "SmallVector move assignment (size 3 to 4)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4 };
      blaze::SmallVector<int,5UL> vec2{ 11, 12, 13 };

      vec2 = std::move( vec1 );

      checkSize    ( vec2, 4UL );
      checkCapacity( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector move assignment (size 8 to 4)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4 };
      blaze::SmallVector<int,5UL> vec2{ 11, 12, 13, 14, 15, 16, 17, 18 };

      vec2 = std::move( vec1 );

      checkSize    ( vec2, 4UL );
      checkCapacity( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector move assignment (size 3 to 5)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5 };
      blaze::SmallVector<int,5UL> vec2{ 11, 12, 13 };

      vec2 = std::move( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector move assignment (size 8 to 5)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5 };
      blaze::SmallVector<int,5UL> vec2{ 11, 12, 13, 14, 15, 16, 17, 18 };

      vec2 = std::move( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector move assignment (size 3 to 6)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5, 6 };
      blaze::SmallVector<int,5UL> vec2{ 11, 12, 13 };

      vec2 = std::move( vec1 );

      checkSize    ( vec2, 6UL );
      checkCapacity( vec2, 6UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 || vec2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector move assignment (size 8 to 6)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5, 6 };
      blaze::SmallVector<int,5UL> vec2{ 11, 12, 13, 14, 15, 16, 17, 18 };

      vec2 = std::move( vec1 );

      checkSize    ( vec2, 6UL );
      checkCapacity( vec2, 6UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 || vec2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SmallVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the SmallVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSubscript()
{
   {
      test_ = "SmallVector::operator[] (size 4)";

      // Assignment to the element at index 2
      blaze::SmallVector<int,5UL> vec{ 0, 0, 1, 0 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      vec[3] = 3;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[2] != 1 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 0
      vec[0] = 4;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 0 1 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 2
      vec[2] += vec[3];

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[2] != 4 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 0 4 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      vec[1] -= 2;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 4 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 3
      vec[3] *= -3;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 4 -9 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 2
      vec[2] /= 2;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 2 || vec[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 2 -9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector::operator[] (size 7)";

      // Assignment to the element at index 2
      blaze::SmallVector<int,5UL> vec{ 0, 0, 1, 0, 0, 0, 0 };

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 5
      vec[5] = 2;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[2] != 1 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      vec[3] = 3;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 0
      vec[0] = 4;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 2
      vec[2] += vec[3];

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[2] != 4 || vec[3] != 3 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 0 4 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      vec[1] -= vec[5];

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != 3 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 4 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 3
      vec[3] *= -3;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != -9 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 4 -9 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 2
      vec[2] /= 2;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 2 || vec[3] != -9 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the SmallVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   {
      test_ = "SmallVector::at() (size 4)";

      // Assignment to the element at index 2
      blaze::SmallVector<int,5UL> vec{ 0, 0, 1, 0 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec.at(2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      vec.at(3) = 3;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[2] != 1 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 0
      vec.at(0) = 4;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 0 1 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 2
      vec.at(2) += vec.at(3);

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[2] != 4 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 0 4 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      vec.at(1) -= 2;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 4 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 3
      vec.at(3) *= -3;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 4 -9 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 2
      vec.at(2) /= 2;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 2 || vec[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 2 -9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector::at() (size 7)";

      // Assignment to the element at index 2
      blaze::SmallVector<int,5UL> vec{ 0, 0, 1, 0, 0, 0, 0 };

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 5
      vec.at(5) = 2;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[2] != 1 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      vec.at(3) = 3;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 0
      vec.at(0) = 4;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 2
      vec.at(2) += vec.at(3);

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[2] != 4 || vec[3] != 3 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 0 4 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      vec.at(1) -= vec.at(5);

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != 3 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 4 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 3
      vec.at(3) *= -3;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != -9 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 4 -9 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 2
      vec.at(2) /= 2;

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 7UL );

      if( vec[0] != 4 || vec[1] != -2 || vec[2] != 2 || vec[3] != -9 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SmallVector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the SmallVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   using VectorType    = blaze::SmallVector<int,5UL>;
   using Iterator      = VectorType::Iterator;
   using ConstIterator = VectorType::ConstIterator;

   VectorType vec{ 1, 0, -2, -3 };

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      Iterator it{};

      if( it != Iterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing the ConstIterator default constructor
   {
      test_ = "ConstIterator default constructor";

      ConstIterator it{};

      if( it != ConstIterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing conversion from Iterator to ConstIterator
   {
      test_ = "Iterator/ConstIterator conversion";

      ConstIterator it( begin( vec ) );

      if( it == end( vec ) || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via Iterator
   {
      test_ = "Iterator subtraction";

      const size_t number( end( vec ) - begin( vec ) );

      if( number != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via ConstIterator
   {
      test_ = "ConstIterator subtraction";

      const size_t number( cend( vec ) - cbegin( vec ) );

      if( number != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      ConstIterator it ( cbegin( vec ) );
      ConstIterator end( cend( vec ) );

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      --it;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it--;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it += 2UL;

      if( it == end || *it != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator addition assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it -= 2UL;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator subtraction assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it + 3UL;

      if( it == end || *it != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar addition failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it - 3UL;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar subtraction failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = 4UL + it;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scalar/iterator addition failed\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing assignment via Iterator
   {
      test_ = "Assignment via Iterator";

      int value = 6;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it = value++;
      }

      if( vec[0] != 6 || vec[1] != 7 || vec[2] != 8 || vec[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      int value = 2;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it += value++;
      }

      if( vec[0] != 8 || vec[1] != 10 || vec[2] != 12 || vec[3] != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 8 10 12 14 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      int value = 2;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it -= value++;
      }

      if( vec[0] != 6 || vec[1] != 7 || vec[2] != 8 || vec[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      int value = 1;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it *= value++;
      }

      if( vec[0] != 6 || vec[1] != 14 || vec[2] != 24 || vec[3] != 36 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 6 14 24 36 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it /= 2;
      }

      if( vec[0] != 3 || vec[1] != 7 || vec[2] != 12 || vec[3] != 18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3 7 12 18 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the SmallVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   test_ = "SmallVector::clear()";

   // Clearing a default constructed vector
   {
      blaze::SmallVector<int,5UL> vec;

      clear( vec );

      checkSize( vec, 0UL );
   }

   // Clearing an initialized vector
   {
      // Initialization check
      blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the vector
      clear( vec );

      checkSize( vec, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the SmallVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   {
      test_ = "SmallVector::resize( size_t )";

      // Initialization check
      blaze::SmallVector<int,5UL> vec;

      checkSize( vec, 0UL );

      // Resizing to 0
      vec.resize( 0UL );

      checkSize( vec, 0UL );

      // Resizing to 4
      vec.resize( 4UL );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 3;
      vec[3] = 4;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6
      vec.resize( 6UL );
      vec[4] = 5;
      vec[5] = 6;

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 3
      vec.resize( 3UL );
      vec[0] = 11;
      vec[1] = 12;
      vec[2] = 13;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );

      if( vec[0] != 11 || vec[1] != 12 || vec[2] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 11 12 13 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6
      vec.resize( 6UL );
      vec[3] = 14;
      vec[4] = 15;
      vec[5] = 16;

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 11 || vec[1] != 12 || vec[2] != 13 || vec[3] != 14 || vec[4] != 15 || vec[5] != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 11 12 13 14 15 16 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0
      vec.resize( 0UL );

      checkSize( vec, 0UL );
   }

   {
      test_ = "SmallVector::resize( size_t, const Type& )";

      // Initialization check
      blaze::SmallVector<int,5UL> vec;

      checkSize( vec, 0UL );

      // Resizing to 0
      vec.resize( 0UL, 2 );

      checkSize( vec, 0UL );

      // Resizing to 4
      vec.resize( 4UL, 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6
      vec.resize( 6UL, 2 );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 || vec[4] != 2 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 3
      vec.resize( 3UL, 2 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6
      vec.resize( 6UL, 2 );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 || vec[4] != 2 || vec[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0
      vec.resize( 0UL, 2 );

      checkSize( vec, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the SmallVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReserve()
{
   test_ = "SmallVector::reserve()";

   // Initialization check
   blaze::SmallVector<int,5UL> vec;

   checkSize( vec, 0UL );

   // Increasing the capacity of the vector
   vec.reserve( 4UL );

   checkSize    ( vec, 0UL );
   checkCapacity( vec, 4UL );

   // Further increasing the capacity of the vector
   vec.reserve( 8UL );

   checkSize    ( vec, 0UL );
   checkCapacity( vec, 8UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member function of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the SmallVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testShrinkToFit()
{
   test_ = "SmallVector::shrinkToFit()";

   // Shrinking a vector without excessive capacity
   {
      blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

      vec.shrinkToFit();

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec.capacity() > 5UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Capacity: " << vec.capacity() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Shrinking a vector with excessive capacity (size 4)
   {
      blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };
      vec.reserve( 100UL );

      vec.shrinkToFit();

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );

      if( vec.capacity() > 5UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Capacity: " << vec.capacity() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Shrinking a vector with excessive capacity (size 8)
   {
      blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6, 7, 8 };
      vec.reserve( 100UL );

      vec.shrinkToFit();

      checkSize    ( vec, 8UL );
      checkCapacity( vec, 8UL );

      if( vec.capacity() > 8UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Capacity: " << vec.capacity() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ||
          vec[4] != 5 || vec[5] != 6 || vec[6] != 7 || vec[7] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c pushBack() member function of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c pushBack() member function of the SmallVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testPushBack()
{
   {
      test_ = "SmallVector::pushBack() (size 4)";

      blaze::SmallVector<int,5UL> vec;

      checkSize( vec, 0UL );

      vec.pushBack( 1 );
      vec.pushBack( 2 );
      vec.pushBack( 3 );
      vec.pushBack( 4 );
      vec.pushBack( 5 );

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      vec.pushBack( 6 );
      vec.pushBack( 7 );
      vec.pushBack( 8 );

      checkSize    ( vec, 8UL );
      checkCapacity( vec, 8UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ||
          vec[4] != 5 || vec[5] != 6 || vec[6] != 7 || vec[7] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the SmallVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testInsert()
{
   {
      // Inserting into an empty small vector
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (empty vector)";

         blaze::SmallVector<int,5UL> vec;
         int value = 1;

         auto pos = vec.insert( vec.begin(), value );

         checkSize    ( vec, 1UL );
         checkCapacity( vec, 1UL );

         if( pos == vec.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small vector (x 2 3 4)
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (x 2 3 4)";

         blaze::SmallVector<int,5UL> vec{ 2, 3, 4 };
         int value = 1;

         auto pos = vec.insert( vec.begin(), value );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos == vec.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small vector (1 x 3 4)
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (1 x 3 4)";

         blaze::SmallVector<int,5UL> vec{ 1, 3, 4 };
         int value = 2;

         auto pos = vec.insert( vec.begin()+1UL, value );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos == vec.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small vector (1 2 3 x)
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (1 2 3 x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3 };
         int value = 4;

         auto pos = vec.insert( vec.end(), value );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos == vec.end() || *pos != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small vector (x 2 3 4 5 6)
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (x 2 3 4 5 6)";

         blaze::SmallVector<int,5UL> vec{ 2, 3, 4, 5, 6 };
         int value = 1;

         auto pos = vec.insert( vec.begin(), value );

         checkSize    ( vec, 6UL );
         checkCapacity( vec, 6UL );

         if( pos == vec.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small vector (1 x 3 4 5 6)
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (1 x 3 4 5 6)";

         blaze::SmallVector<int,5UL> vec{ 1, 3, 4, 5, 6 };
         int value = 2;

         auto pos = vec.insert( vec.begin()+1UL, value );

         checkSize    ( vec, 6UL );
         checkCapacity( vec, 6UL );

         if( pos == vec.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small vector (1 2 3 4 5 x)
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (1 2 3 4 5 x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5 };
         int value = 6;

         auto pos = vec.insert( vec.end(), value );

         checkSize    ( vec, 6UL );
         checkCapacity( vec, 6UL );

         if( pos == vec.end() || *pos != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 6\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small vector (x 2 3 4 5 6 7 8)
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (x 2 3 4 5 6 7 8)";

         blaze::SmallVector<int,5UL> vec{ 2, 3, 4, 5, 6, 7, 8 };
         int value = 1;

         auto pos = vec.insert( vec.begin(), value );

         checkSize    ( vec, 8UL );
         checkCapacity( vec, 8UL );

         if( pos == vec.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ||
             vec[4] != 5 || vec[5] != 6 || vec[6] != 7 || vec[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small vector (1 x 3 4 5 6 7 8 )
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (1 x 3 4 5 6 7 8)";

         blaze::SmallVector<int,5UL> vec{ 1, 3, 4, 5, 6, 7, 8 };
         int value = 2;

         auto pos = vec.insert( vec.begin()+1UL, value );

         checkSize    ( vec, 8UL );
         checkCapacity( vec, 8UL );

         if( pos == vec.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ||
             vec[4] != 5 || vec[5] != 6 || vec[6] != 7 || vec[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small vector (1 2 3 4 5 6 7 x)
      {
         test_ = "SmallVector::insert( Iterator, const Type& ) (1 2 3 4 5 6 7 x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6, 7 };
         int value = 8;

         auto pos = vec.insert( vec.end(), value );

         checkSize    ( vec, 8UL );
         checkCapacity( vec, 8UL );

         if( pos == vec.end() || *pos != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 8\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ||
             vec[4] != 5 || vec[5] != 6 || vec[6] != 7 || vec[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      // Inserting into an empty small vector
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (empty vector)";

         blaze::SmallVector<int,5UL> vec;

         auto pos = vec.insert( vec.begin(), 1 );

         checkSize    ( vec, 1UL );
         checkCapacity( vec, 1UL );

         if( pos == vec.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small vector (x 2 3 4)
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (x 2 3 4)";

         blaze::SmallVector<int,5UL> vec{ 2, 3, 4 };

         auto pos = vec.insert( vec.begin(), 1 );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos == vec.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small vector (1 x 3 4)
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (1 x 3 4)";

         blaze::SmallVector<int,5UL> vec{ 1, 3, 4 };

         auto pos = vec.insert( vec.begin()+1UL, 2 );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos == vec.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small vector (1 2 3 x)
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (1 2 3 x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3 };

         auto pos = vec.insert( vec.end(), 4 );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos == vec.end() || *pos != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small vector (x 2 3 4 5 6)
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (x 2 3 4 5 6)";

         blaze::SmallVector<int,5UL> vec{ 2, 3, 4, 5, 6 };

         auto pos = vec.insert( vec.begin(), 1 );

         checkSize    ( vec, 6UL );
         checkCapacity( vec, 6UL );

         if( pos == vec.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small vector (1 x 3 4 5 6)
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (1 x 3 4 5 6)";

         blaze::SmallVector<int,5UL> vec{ 1, 3, 4, 5, 6 };

         auto pos = vec.insert( vec.begin()+1UL, 2 );

         checkSize    ( vec, 6UL );
         checkCapacity( vec, 6UL );

         if( pos == vec.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small vector (1 2 3 4 5 x)
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (1 2 3 4 5 x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5 };

         auto pos = vec.insert( vec.end(), 6 );

         checkSize    ( vec, 6UL );
         checkCapacity( vec, 6UL );

         if( pos == vec.end() || *pos != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 6\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 || vec[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small vector (x 2 3 4 5 6 7 8)
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (x 2 3 4 5 6 7 8)";

         blaze::SmallVector<int,5UL> vec{ 2, 3, 4, 5, 6, 7, 8 };

         auto pos = vec.insert( vec.begin(), 1 );

         checkSize    ( vec, 8UL );
         checkCapacity( vec, 8UL );

         if( pos == vec.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ||
             vec[4] != 5 || vec[5] != 6 || vec[6] != 7 || vec[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small vector (1 x 3 4 5 6 7 8 )
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (1 x 3 4 5 6 7 8)";

         blaze::SmallVector<int,5UL> vec{ 1, 3, 4, 5, 6, 7, 8 };

         auto pos = vec.insert( vec.begin()+1UL, 2 );

         checkSize    ( vec, 8UL );
         checkCapacity( vec, 8UL );

         if( pos == vec.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ||
             vec[4] != 5 || vec[5] != 6 || vec[6] != 7 || vec[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small vector (1 2 3 4 5 6 7 x)
      {
         test_ = "SmallVector::insert( Iterator, Type&& ) (1 2 3 4 5 6 7 x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6, 7 };

         auto pos = vec.insert( vec.end(), 8 );

         checkSize    ( vec, 8UL );
         checkCapacity( vec, 8UL );

         if( pos == vec.end() || *pos != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 8\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ||
             vec[4] != 5 || vec[5] != 6 || vec[6] != 7 || vec[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the SmallVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testErase()
{
   {
      // Erasing from the beginning of a small vector
      {
         test_ = "SmallVector::erase( Iterator ) (x 2 3 4)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

         auto pos = vec.erase( vec.begin() );

         checkSize    ( vec, 3UL );
         checkCapacity( vec, 3UL );

         if( pos == vec.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 2 || vec[1] != 3 || vec[2] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the middle of a small vector
      {
         test_ = "SmallVector::erase( Iterator ) (1 x 3 4)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

         auto pos = vec.erase( vec.begin()+1UL );

         checkSize    ( vec, 3UL );
         checkCapacity( vec, 3UL );

         if( pos == vec.end() || *pos != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 3 || vec[2] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the end of a small vector
      {
         test_ = "SmallVector::erase( Iterator ) (1 2 3 x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

         auto pos = vec.erase( vec.begin()+3UL );

         checkSize    ( vec, 3UL );
         checkCapacity( vec, 3UL );

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the beginning of a small vector
      {
         test_ = "SmallVector::erase( Iterator ) (x 2 3 4 5 6)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6 };

         auto pos = vec.erase( vec.begin() );

         checkSize    ( vec, 5UL );
         checkCapacity( vec, 5UL );

         if( pos == vec.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 2 || vec[1] != 3 || vec[2] != 4 || vec[3] != 5 || vec[4] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the middle of a small vector
      {
         test_ = "SmallVector::erase( Iterator ) (1 2 x 4 5 6)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6 };

         auto pos = vec.erase( vec.begin()+2UL );

         checkSize    ( vec, 5UL );
         checkCapacity( vec, 5UL );

         if( pos == vec.end() || *pos != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 4 || vec[3] != 5 || vec[4] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the end of a small vector
      {
         test_ = "SmallVector::erase( Iterator ) (1 2 3 4 5 x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6 };

         auto pos = vec.erase( vec.begin()+5UL );

         checkSize    ( vec, 5UL );
         checkCapacity( vec, 5UL );

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 || vec[4] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 5 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      // Erasing from the beginning of a small vector
      {
         test_ = "SmallVector::erase( Iterator, Iterator ) (x x 3 4)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

         auto pos = vec.erase( vec.begin(), vec.begin()+2UL );

         checkSize    ( vec, 2UL );
         checkCapacity( vec, 2UL );

         if( pos == vec.end() || *pos != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 3 || vec[1] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the middle of a small vector
      {
         test_ = "SmallVector::erase( Iterator, Iterator ) (1 x x 4)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

         auto pos = vec.erase( vec.begin()+1UL, vec.begin()+3UL );

         checkSize    ( vec, 2UL );
         checkCapacity( vec, 2UL );

         if( pos == vec.end() || *pos != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the end of a small vector
      {
         test_ = "SmallVector::erase( Iterator, Iterator ) (1 2 x x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4 };

         auto pos = vec.erase( vec.begin()+2UL, vec.begin()+4UL );

         checkSize    ( vec, 2UL );
         checkCapacity( vec, 2UL );

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the beginning of a small vector
      {
         test_ = "SmallVector::erase( Iterator, Iterator ) (x x 3 4 5 6)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6 };

         auto pos = vec.erase( vec.begin(), vec.begin()+2UL );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos == vec.end() || *pos != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 3 || vec[1] != 4 || vec[2] != 5 || vec[3] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the middle of a small vector
      {
         test_ = "SmallVector::erase( Iterator, Iterator ) (1 2 x x 5 6)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6 };

         auto pos = vec.erase( vec.begin()+2UL, vec.begin()+4UL );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos == vec.end() || *pos != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 5\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 5 || vec[3] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the end of a small vector
      {
         test_ = "SmallVector::erase( Iterator, Iterator ) (1 2 3 4 x x)";

         blaze::SmallVector<int,5UL> vec{ 1, 2, 3, 4, 5, 6 };

         auto pos = vec.erase( vec.begin()+4UL, vec.begin()+6UL );

         checkSize    ( vec, 4UL );
         checkCapacity( vec, 4UL );

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}

//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the SmallVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the SmallVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   {
      test_ = "SmallVector swap (size 3 and size 4)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3 };
      blaze::SmallVector<int,5UL> vec2{ 4, 3, 2, 1 };

      swap( vec1, vec2 );

      checkSize    ( vec1, 4UL );
      checkCapacity( vec1, 4UL );

      if( vec1[0] != 4 || vec1[1] != 3 || vec1[2] != 2 || vec1[3] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 4 3 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkSize    ( vec2, 3UL );
      checkCapacity( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector swap (size 3 and size 7)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3 };
      blaze::SmallVector<int,5UL> vec2{ 7, 6, 5, 4, 3, 2, 1 };

      swap( vec1, vec2 );

      checkSize    ( vec1, 7UL );
      checkCapacity( vec1, 7UL );

      if( vec1[0] != 7 || vec1[1] != 6 || vec1[2] != 5 || vec1[3] != 4 ||
          vec1[4] != 3 || vec1[5] != 2 || vec1[6] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 7 6 5 4 3 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkSize    ( vec2, 3UL );
      checkCapacity( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector swap (size 8 and size 4)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5, 6, 7, 8 };
      blaze::SmallVector<int,5UL> vec2{ 4, 3, 2, 1 };

      swap( vec1, vec2 );

      checkSize    ( vec1, 4UL );
      checkCapacity( vec1, 4UL );

      if( vec1[0] != 4 || vec1[1] != 3 || vec1[2] != 2 || vec1[3] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 4 3 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkSize    ( vec2, 8UL );
      checkCapacity( vec2, 8UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ||
          vec2[4] != 5 || vec2[5] != 6 || vec2[6] != 7 || vec2[7] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallVector swap (size 8 and size 7)";

      blaze::SmallVector<int,5UL> vec1{ 1, 2, 3, 4, 5, 6, 7, 8 };
      blaze::SmallVector<int,5UL> vec2{ 7, 6, 5, 4, 3, 2, 1 };

      swap( vec1, vec2 );

      checkSize    ( vec1, 7UL );
      checkCapacity( vec1, 7UL );

      if( vec1[0] != 7 || vec1[1] != 6 || vec1[2] != 5 || vec1[3] != 4 ||
          vec1[4] != 3 || vec1[5] != 2 || vec1[6] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 7 6 5 4 3 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkSize    ( vec2, 8UL );
      checkCapacity( vec2, 8UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ||
          vec2[4] != 5 || vec2[5] != 6 || vec2[6] != 7 || vec2[7] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace smallvector

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
   std::cout << "   Running SmallVector class test..." << std::endl;

   try
   {
      RUN_SMALLVECTOR_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SmallVector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
