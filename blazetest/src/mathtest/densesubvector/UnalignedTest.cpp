//=================================================================================================
/*!
//  \file src/mathtest/densesubvector/UnalignedTest.cpp
//  \brief Source file for the unaligned DenseSubvector class test
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <blaze/math/CompressedVector.h>
#include <blazetest/mathtest/densesubvector/UnalignedTest.h>


namespace blazetest {

namespace mathtest {

namespace densesubvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the unaligned DenseSubvector class test.
//
// \exception std::runtime_error Operation error detected.
*/
UnalignedTest::UnalignedTest()
   : vec_( 8UL )
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testDivAssign();
   testSubscript();
   testIterator();
   testNonZeros();
   testReset();
   testScale();
   testIsDefault();
   testIsNan();
   testMinimum();
   testMaximum();
   testSubvector();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the DenseSubvector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the DenseSubvector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testConstructors()
{
   test_ = "DenseSubvector constructor";

   initialize();

   for( size_t start=0UL; start<vec_.size(); ++start ) {
      for( size_t size=1UL; start+size<vec_.size(); ++size )
      {
         SVT sv = subvector( vec_, start, size );

         for( size_t i=0UL; i<size; ++i )
         {
            if( sv[i] != vec_[start+i] ) {
               std::ostringstream oss;
               oss << " Test: " << test_ << "\n"
                   << " Error: Setup of dense subvector failed\n"
                   << " Details:\n"
                   << "   Start = " << start << "\n"
                   << "   Size  = " << size  << "\n"
                   << "   Subvector:\n" << sv << "\n"
                   << "   Vector:\n" << vec_ << "\n";
               throw std::runtime_error( oss.str() );
            }
         }
      }
   }

   try {
      SVT sv = subvector( vec_, 2UL, 7UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      SVT sv = subvector( vec_, 9UL, 0UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DenseSubvector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the DenseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testAssignment()
{
   //=====================================================================================
   // Homogeneous assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector homogeneous assignment";

      initialize();

      SVT sv = subvector( vec_, 2UL, 4UL );
      sv = 12;

      checkSize    ( sv  ,  4UL );
      checkNonZeros( sv  ,  4UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  6UL );

      if( sv[0] != 12 || sv[1] != 12 || sv[2] != 12 || sv[3] != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 12 12 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] !=  1 || vec_[2] != 12 || vec_[3] != 12 ||
          vec_[4] != 12 || vec_[5] != 12 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 12 12 12 12 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector copy assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      SVT sv = subvector( vec, 5UL, 3UL );
      sv = subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  2UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  2UL );

      if( sv[0] != -3 || sv[1] != 0 || sv[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] !=  0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != -3 || vec[6] != 0 || vec[7] != 4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 -3 0 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "DenseSubvector copy assignment (aliasing)";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );
      sv = subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != -3 || sv[1] != 0 || sv[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -3 || vec_[2] != 0 || vec_[3] != 4 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -3 0 4 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "Dense vector assignment";

      initialize();

      SVT sv = subvector( vec_, 3UL, 4UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );
      vec[1] = 8;
      vec[3] = 9;

      sv = vec;

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv != vec ||
          sv[0] != 0 || sv[1] != 8 || sv[2] != 0 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != 8 || vec_[5] != 0 || vec_[6] != 9 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 0 8 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "Sparse vector assignment";

      initialize();

      SVT sv = subvector( vec_, 3UL, 4UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL, 1UL );
      vec[3] = 9;

      sv = vec;

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv != vec ||
          sv[0] != 0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != 0 || vec_[5] != 0 || vec_[6] != 9 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 0 0 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DenseSubvector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the DenseSubvector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testAddAssign()
{
   //=====================================================================================
   // DenseSubvector addition assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector addition assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      SVT sv = subvector( vec, 5UL, 3UL );
      sv += subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  3UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  3UL );

      if( sv[0] != 3 || sv[1] != -8 || sv[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 -8 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 0 || vec[1] !=  0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 3 || vec[6] != -8 || vec[7] != 4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 0 0 3 -8 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "DenseSubvector addition assignment (aliasing)";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );
      sv += subvector( vec_, 3UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -1 || sv[1] != -3 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -1 -3 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -1 || vec_[2] != -3 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -1 -3 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Dense vector addition assignment";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      sv += vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != 3 || sv[1] != -4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 -4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != -4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 -4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Sparse vector addition assignment";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] =  2;
      vec[1] = -4;

      sv += vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != 3 || sv[1] != -4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 -4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != -4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 -4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DenseSubvector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the DenseSubvector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testSubAssign()
{
   //=====================================================================================
   // DenseSubvector subtraction assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector subtraction assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      SVT sv = subvector( vec, 5UL, 3UL );
      sv -= subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  3UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  3UL );

      if( sv[0] != 9 || sv[1] != -8 || sv[2] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 -8 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 0 || vec[1] !=  0 || vec[2] !=  0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 9 || vec[6] != -8 || vec[7] != -4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 0 0 9 -8 -4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "DenseSubvector subtraction assignment (aliasing)";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );
      sv -= subvector( vec_, 3UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != 3 || sv[1] != 3 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 3 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 3 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 3 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Dense vector subtraction assignment";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      sv -= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -1 || sv[1] != 4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -1 4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -1 || vec_[2] != 4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -1 4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Sparse vector subtraction assignment";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] =  2;
      vec[1] = -4;

      sv -= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -1 || sv[1] != 4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -1 4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -1 || vec_[2] != 4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -1 4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DenseSubvector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the DenseSubvector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testMultAssign()
{
   //=====================================================================================
   // DenseSubvector multiplication assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector multiplication assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      SVT sv = subvector( vec, 5UL, 3UL );
      sv *= subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  1UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  1UL );

      if( sv[0] != -18 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -18 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] !=   0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != -18 || vec[6] != 0 || vec[7] != 0 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 0 0 -18 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "DenseSubvector multiplication assignment (aliasing)";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );
      sv *= subvector( vec_, 3UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != -2 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Dense vector multiplication assignment";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      sv *= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Sparse vector multiplication assignment";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] =  2;
      vec[1] = -4;

      sv *= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Scalar multiplication assignment";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );

      sv *= 3;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 3 || sv[1] != 0 || sv[2] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DenseSubvector division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the DenseSubvector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testDivAssign()
{
   //=====================================================================================
   // Scalar division assignment
   //=====================================================================================

   {
      test_ = "Scalar division assignment";

      initialize();

      SVT sv = subvector( vec_, 1UL, 3UL );

      sv /= 0.5;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != -4 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 -4 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DenseSubvector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the DenseSubvector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void UnalignedTest::testSubscript()
{
   test_ = "DenseSubvector::operator[]";

   initialize();

   SVT sv = subvector( vec_, 1UL, 4UL );

   // Writing the first element
   sv[1] = 9;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 4UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 5UL );

   if( sv[0] != 1 || sv[1] != 9 || sv[2] != -2 || sv[3] != -3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 1 9 -2 -3 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != -2 ||
       vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 -2 -3 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Writing the second element
   sv[2] = 0;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 1 || sv[1] != 9 || sv[2] != 0 || sv[3] != -3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 1 9 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != 0 ||
       vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 0 -3 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Writing the third element
   sv[3] = -8;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 1 || sv[1] != 9 || sv[2] != 0 || sv[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( -2 9 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DenseSubvector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the DenseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testIterator()
{
   initialize();

   // Counting the number of elements in first half of the vector
   {
      test_ = "Iterator subtraction";

      SVT sv = subvector( vec_, 0UL, 5UL );
      const size_t number( sv.end() - sv.begin() );

      if( number != 5UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 5\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in second half of the vector
   {
      test_ = "Iterator subtraction";

      SVT sv = subvector( vec_, 5UL, 3UL );
      const size_t number( sv.end() - sv.begin() );

      if( number != 3UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 3\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      SVT sv = subvector( vec_, 1UL, 4UL );
      SVT::ConstIterator it ( sv.cbegin() );
      SVT::ConstIterator end( sv.cend() );

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

      SVT sv = subvector( vec_, 2UL, 4UL );
      int value = 6;

      for( SVT::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
         *it = value++;
      }

      if( sv[0] != 6 || sv[1] != 7 || sv[2] != 8 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 6 || vec_[3] != 7 ||
          vec_[4] != 8 || vec_[5] != 9 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 6 7 8 9 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      SVT sv = subvector( vec_, 2UL, 4UL );
      int value = 2;

      for( SVT::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
         *it += value++;
      }

      if( sv[0] != 8 || sv[1] != 10 || sv[2] != 12 || sv[3] != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 8 10 12 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] !=  1 || vec_[2] != 8 || vec_[3] != 10 ||
          vec_[4] != 12 || vec_[5] != 14 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 8 10 12 14 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      SVT sv = subvector( vec_, 2UL, 4UL );
      int value = 2;

      for( SVT::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
         *it -= value++;
      }

      if( sv[0] != 6 || sv[1] != 7 || sv[2] != 8 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 6 || vec_[3] != 7 ||
          vec_[4] != 8 || vec_[5] != 9 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 6 7 8 9 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      SVT sv = subvector( vec_, 2UL, 4UL );
      int value = 1;

      for( SVT::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
         *it *= value++;
      }

      if( sv[0] != 6 || sv[1] != 14 || sv[2] != 24 || sv[3] != 36 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 6 14 24 36 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] !=  1 || vec_[2] != 6 || vec_[3] != 14 ||
          vec_[4] != 24 || vec_[5] != 36 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 6 14 24 36 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      SVT sv = subvector( vec_, 2UL, 4UL );

      for( SVT::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
         *it /= 2;
      }

      if( sv[0] != 3 || sv[1] != 7 || sv[2] != 12 || sv[3] != 18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 7 12 18 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] !=  1 || vec_[2] != 3 || vec_[3] != 7 ||
          vec_[4] != 12 || vec_[5] != 18 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 3 7 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the nonZeros member function of DenseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of DenseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testNonZeros()
{
   test_ = "DenseSubvector::nonZeros()";

   initialize();

   // Initialization check
   SVT sv = subvector( vec_, 0UL, 4UL );

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 2UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 1 0 -2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense subvector
   sv[3] = 0;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 1UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 3UL );

   if( sv[0] != 0 || sv[1] != 1 || sv[2] != 0 || sv[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 1 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense vector
   vec_[2UL] = 5;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 2UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 0 || sv[1] != 1 || sv[2] != 5 || sv[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 1 5 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reset member function of DenseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of DenseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testReset()
{
   test_ = "DenseSubvector::reset()";

   initialize();

   // Resetting the range [0,3]
   {
      SVT sv = subvector( vec_, 0UL, 4UL );
      sv.reset();

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 0UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv[0] != 0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting the range [4,7]
   {
      SVT sv = subvector( vec_, 4UL, 4UL );
      sv.reset();

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 0UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 0UL );

      if( sv[0] != 0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [4,7] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of DenseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of DenseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testScale()
{
   test_ = "DenseSubvector::scale()";

   initialize();

   SVT sv = subvector( vec_, 1UL, 4UL );

   // Integral scaling of the subvector in the range [1,4]
   sv.scale( 3 );

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 3 || sv[1] != 0 || sv[2] != -6 || sv[3] != -9 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Integral scale operation of range [1,4] failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 3 0 -6 -9 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
       vec_[4] != -9 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Integral scale operation of range [1,4] failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 3 0 -6 -9 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Floating point scaling of the subvector in the range [1,4]
   sv.scale( 0.5 );

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 1 || sv[1] != 0 || sv[2] != -3 || sv[3] != -4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Floating point scale operation of range [1,4] failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 1 0 -3 -4 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -3 ||
       vec_[4] != -4 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Floating point scale operation of range [1,4] failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 0 -3 -4 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDefault function with the DenseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDefault function with the DenseSubvector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testIsDefault()
{
   test_ = "isDefault() function";

   initialize();

   // isDefault with default vector
   {
      VT vec( 8UL, 0 );
      SVT sv = subvector( vec, 2UL, 5UL );

      if( isDefault( sv ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with non-default vector
   {
      SVT sv = subvector( vec_, 2UL, 5UL );

      if( isDefault( sv ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isnan function with the DenseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isnan function with the DenseSubvector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testIsNan()
{
   test_ = "isnan() function";

   typedef blaze::DynamicVector<float,blaze::columnVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType>                SubvectorType;

   VectorType vec( 9UL, 0 );
   vec[2] =  1;
   vec[3] = -2;
   vec[4] = -3;
   vec[8] =  4;

   // isnan with empty 3-dimensional subvector
   {
      SubvectorType sv = subvector( vec, 5UL, 3UL );

      checkSize    ( sv, 3UL );
      checkNonZeros( sv, 0UL );

      if( blaze::isnan( sv ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with partially filled 5-dimensional subvector
   {
      SubvectorType sv = subvector( vec, 4UL, 5UL );

      checkSize    ( sv, 5UL );
      checkNonZeros( sv, 2UL );

      if( blaze::isnan( sv ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with fully filled 3-dimensional subvector
   {
      SubvectorType sv = subvector( vec, 2UL, 3UL );

      checkSize    ( sv, 3UL );
      checkNonZeros( sv, 3UL );

      if( blaze::isnan( sv ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the min function with the DenseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the min function used with the DenseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testMinimum()
{
   test_ = "min() function";

   initialize();

   // Computing the minimum of the in the range [0,2]
   {
      const int minimum = min( subvector( vec_, 0UL, 3UL ) );

      if( minimum != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Minimum computation for range [0,2] failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the minimum of the in the range [2,4]
   {
      const int minimum = min( subvector( vec_, 2UL, 3UL ) );

      if( minimum != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Minimum computation for range [2,4] failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the minimum of the in the range [4,6]
   {
      const int minimum = min( subvector( vec_, 4UL, 3UL ) );

      if( minimum != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Minimum computation for range [4,6] failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the minimum of the in the range [6,7]
   {
      const int minimum = min( subvector( vec_, 6UL, 2UL ) );

      if( minimum != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Minimum computation for range [6,7] failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the max function with the DenseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the max function used with the DenseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testMaximum()
{
   test_ = "max() function";

   initialize();

   // Computing the maximum of the in the range [0,2]
   {
      const int maximum = max( subvector( vec_, 0UL, 3UL ) );

      if( maximum != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Maximum computation for range [0,2] failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the maximum of the in the range [2,4]
   {
      const int maximum = max( subvector( vec_, 2UL, 3UL ) );

      if( maximum != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Maximum computation for range [2,4] failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the maximum of the in the range [4,6]
   {
      const int maximum = max( subvector( vec_, 4UL, 3UL ) );

      if( maximum != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Maximum computation for range [4,6] failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the maximum of the in the range [6,7]
   {
      const int maximum = max( subvector( vec_, 6UL, 2UL ) );

      if( maximum != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Maximum computation for range [6,7] failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subvector function with the DenseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subvector function used with the DenseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void UnalignedTest::testSubvector()
{
   test_ = "subvector() function";

   initialize();

   {
      SVT sv1 = subvector( vec_, 1UL, 6UL );
      SVT sv2 = subvector( sv1 , 1UL, 4UL );

      if( sv2[1] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << sv2[1] << "\n"
             << "   Expected result: -2\n";
         throw std::runtime_error( oss.str() );
      }

      if( *sv2.begin() != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *sv2.begin() << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      SVT sv1 = subvector( vec_, 1UL, 6UL );
      SVT sv2 = subvector( sv1 , 6UL, 2UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      SVT sv1 = subvector( vec_, 1UL, 6UL );
      SVT sv2 = subvector( sv1 , 2UL, 5UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member vectors to specific predetermined values.
*/
void UnalignedTest::initialize()
{
   // Initializing the dynamic row vector
   vec_[0] =  0;
   vec_[1] =  1;
   vec_[2] =  0;
   vec_[3] = -2;
   vec_[4] = -3;
   vec_[5] =  0;
   vec_[6] =  4;
   vec_[7] =  0;
}
//*************************************************************************************************

} // namespace densesubvector

} // namespace mathtest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running unaligned DenseSubvector class test..." << std::endl;

   try
   {
      RUN_DENSESUBVECTOR_UNALIGNED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during unaligned DenseSubvector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
