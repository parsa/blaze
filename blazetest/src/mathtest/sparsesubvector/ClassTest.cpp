//=================================================================================================
/*!
//  \file src/mathtest/sparsesubvector/ClassTest.cpp
//  \brief Source file for the SparseSubvector class test
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
#include <blaze/math/DynamicVector.h>
#include <blazetest/mathtest/sparsesubvector/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace sparsesubvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SparseSubvector class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
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
   testAppend();
   testInsert();
   testErase();
   testScale();
   testFind();
   testLowerBound();
   testUpperBound();
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
/*!\brief Test of the SparseSubvector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the SparseSubvector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   test_ = "SparseSubvector constructor";

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
                   << " Error: Setup of sparse subvector failed\n"
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
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseSubvector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the SparseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "SparseSubvector copy assignment (no aliasing)";

      initialize();

      VT vec( 10UL );
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
      test_ = "SparseSubvector copy assignment (aliasing)";

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL );
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
/*!\brief Test of the SparseSubvector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the SparseSubvector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // SparseSubvector addition assignment
   //=====================================================================================

   {
      test_ = "SparseSubvector addition assignment (no aliasing)";

      initialize();

      VT vec( 10UL );
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
      test_ = "SparseSubvector addition assignment (aliasing)";

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL );
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
/*!\brief Test of the SparseSubvector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the SparseSubvector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // SparseSubvector subtraction assignment
   //=====================================================================================

   {
      test_ = "SparseSubvector subtraction assignment (no aliasing)";

      initialize();

      VT vec( 10UL );
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
      test_ = "SparseSubvector subtraction assignment (aliasing)";

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL );
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
/*!\brief Test of the SparseSubvector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the SparseSubvector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // SparseSubvector multiplication assignment
   //=====================================================================================

   {
      test_ = "SparseSubvector multiplication assignment (no aliasing)";

      initialize();

      VT vec( 10UL );
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
      test_ = "SparseSubvector multiplication assignment (aliasing)";

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 0 );
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
/*!\brief Test of the SparseSubvector division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the SparseSubvector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testDivAssign()
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
/*!\brief Test of the SparseSubvector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the SparseSubvector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSubscript()
{
   test_ = "SparseSubvector::operator[]";

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
/*!\brief Test of the SparseSubvector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the SparseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   initialize();

   // Counting the number of elements in first half of the vector
   {
      test_ = "Iterator subtraction";

      SVT sv = subvector( vec_, 0UL, 4UL );
      const size_t number( sv.end() - sv.begin() );

      if( number != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in second half of the vector
   {
      test_ = "Iterator subtraction";

      SVT sv = subvector( vec_, 4UL, 4UL );
      const size_t number( sv.end() - sv.begin() );

      if( number != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      SVT sv = subvector( vec_, 1UL, 3UL );
      SVT::ConstIterator it ( sv.cbegin() );
      SVT::ConstIterator end( sv.cend() );

      if( it == end || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || it->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
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

      if( sv[0] != 0 || sv[1] != 6 || sv[2] != 7 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 6 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 7 0 4 0 )\n";
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

      if( sv[0] != 0 || sv[1] != 8 || sv[2] != 10 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 10 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 8 ||
          vec_[4] != 10 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 8 10 0 4 0 )\n";
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

      if( sv[0] != 0 || sv[1] != 6 || sv[2] != 7 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 10 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 7 0 4 0 )\n";
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

      if( sv[0] != 0 || sv[1] != 6 || sv[2] != 14 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 10 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 14 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 14 0 4 0 )\n";
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

      if( sv[0] != 0 || sv[1] != 3 || sv[2] != 7 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 10 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 3 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
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
/*!\brief Test of the nonZeros member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of SparseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   test_ = "SparseSubvector::nonZeros()";

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

   // Changing the number of non-zeros via the sparse subvector
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

   // Changing the number of non-zeros via the sparse vector
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
/*!\brief Test of the reset member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of SparseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   test_ = "SparseSubvector::reset()";

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
/*!\brief Test of the append member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the append member function of SparseSubvector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAppend()
{
   test_ = "SparseSubvector::append()";

   VT vec( 10UL );

   SVT sv = subvector( vec, 2UL, 4UL );
   sv.reserve( 4UL );

   // Appending one non-zero element
   sv.append( 0UL, 1 );

   checkSize    ( sv , 4UL );
   checkCapacity( sv , 4UL );
   checkNonZeros( sv , 1UL );
   checkNonZeros( vec, 1UL );

   if( sv[0] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 0 1 0 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Appending three more non-zero elements
   sv.append( 1UL, 2 );
   sv.append( 2UL, 3 );
   sv.append( 3UL, 4 );

   checkSize    ( sv , 4UL );
   checkCapacity( sv , 4UL );
   checkNonZeros( sv , 4UL );
   checkNonZeros( vec, 4UL );

   if( sv[0] != 1 || sv[1] != 2 || sv[2] != 3 || sv[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Append operation failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 0 1 2 3 4 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the insert member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the insert member function of SparseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testInsert()
{
   test_ = "SparseSubvector::insert()";

   initialize();

   SVT sv = subvector( vec_, 0UL, 8UL );

   // Inserting a non-zero element at the end of the subvector
   {
      SVT::Iterator pos = sv.insert( 7UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 5UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( pos->value() != 9 || pos->index() != 7UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 7\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a non-zero element at the beginning of the subvector
   {
      SVT::Iterator pos = sv.insert( 0UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 6UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 6UL );

      if( pos->value() != 9 || pos->index() != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 0\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  9 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a non-zero element at the center of the subvector
   {
      SVT::Iterator pos = sv.insert( 2UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 7UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 7UL );

      if( pos->value() != 9 || pos->index() != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 2\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  9 || sv[1] != 1 || sv[2] != 9 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 2 9 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Trying to insert an already existing element
   try {
      sv.insert( 3UL, 9 );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Inserting an existing element succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 9 2 0 9 -3 0 4 9 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the erase member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the erase member function of SparseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testErase()
{
   //=====================================================================================
   // Index-based erase function
   //=====================================================================================

   {
      test_ = "SparseSubvector::erase( size_t )";

      initialize();

      SVT sv = subvector( vec_, 1UL, 6UL );

      // Erasing the non-zero element at the end of the subvector
      sv.erase( 5UL );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] !=  1 || sv[1] != 0 || sv[2] != -2 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the subvector
      sv.erase( 0UL );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != -2 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the subvector
      sv.erase( 2UL );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 1UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      sv.erase( 2UL );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 1UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Iterator-based erase function
   //=====================================================================================

   {
      test_ = "SparseSubvector::erase( Iterator )";

      initialize();

      SVT sv = subvector( vec_, 1UL, 6UL );

      // Erasing the non-zero element at the end of the subvector
      {
         SVT::Iterator pos = sv.erase( sv.find( 5UL ) );

         checkSize    ( sv  , 6UL );
         checkNonZeros( sv  , 3UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 3UL );

         if( pos != sv.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  1 || sv[1] != 0 || sv[2] != -2 ||
             sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 1 0 -2 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the subvector
      {
         SVT::Iterator pos = sv.erase( sv.find( 0UL ) );

         checkSize    ( sv  , 6UL );
         checkNonZeros( sv  , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos->value() != -2 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -2\n"
                << "   Expected index:  2\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 0 || sv[2] != -2 ||
             sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the subvector
      {
         SVT::Iterator pos = sv.erase( sv.find( 2UL ) );

         checkSize    ( sv  , 6UL );
         checkNonZeros( sv  , 1UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 1UL );

         if( pos->value() != -3 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -3\n"
                << "   Expected index:  3\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 ||
             sv[3] != -3 || sv[4] != 0 || sv[5] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         SVT::Iterator pos = sv.erase( sv.find( 2UL ) );

         checkSize    ( sv  , 6UL );
         checkNonZeros( sv  , 1UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 1UL );

         if( pos != sv.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 ||
             sv[3] != -3 || sv[4] != 0 || sv[5] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "SparseSubvector::erase( Iterator, Iterator )";

      // Erasing the entire vector
      {
         initialize();

         SVT sv = subvector( vec_, 0UL, 8UL );

         SVT::Iterator pos = sv.erase( sv.begin(), sv.end() );

         checkSize    ( sv  , 8UL );
         checkNonZeros( sv  , 0UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 0UL );

         if( pos != sv.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] != 0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 0 ||
             sv[4] != 0 || sv[5] != 0 || sv[6] != 0 || sv[7] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the subvector failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 0 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the vector
      {
         initialize();

         SVT sv = subvector( vec_, 0UL, 8UL );

         SVT::Iterator pos = sv.erase( sv.begin(), sv.find( 4UL ) );

         checkSize    ( sv  , 8UL );
         checkNonZeros( sv  , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos->value() != -3 || pos->index() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -3\n"
                << "   Expected index:  4\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 0 ||
             sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial subvector failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 0 0 -3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the vector
      {
         initialize();

         SVT sv = subvector( vec_, 0UL, 8UL );

         SVT::Iterator pos = sv.erase( sv.find( 4UL ), sv.end() );

         checkSize    ( sv  , 8UL );
         checkNonZeros( sv  , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos != sv.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] != 0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
             sv[4] != 0 || sv[5] != 0 || sv[6] != 0 || sv[7] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial subvector failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 1 0 -2 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         initialize();

         SVT sv = subvector( vec_, 0UL, 8UL );

         SVT::Iterator pos = sv.erase( sv.find( 1UL ), sv.find( 1UL ) );

         checkSize    ( sv  , 8UL );
         checkNonZeros( sv  , 4UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 4UL );

         if( pos != sv.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
             sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 1 0 -2 -3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of SparseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScale()
{
   test_ = "SparseSubvector::scale()";

   initialize();

   SVT sv = subvector( vec_, 1UL, 4UL );

   // Integral scaling the subvector in the range [1,4]
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

   // Floating point scaling the subvector in the range [1,4]
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
/*!\brief Test of the find member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the find member function of SparseSubvector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   test_ = "SparseSubvector::find()";

   initialize();

   SVT sv = subvector( vec_, 1UL, 5UL );

   // Searching for the first element
   {
      SVT::Iterator pos = sv.find( 0UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 0 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the second element
   {
      SVT::Iterator pos = sv.find( 2UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = -2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for a non-existing non-zero element
   {
      SVT::Iterator pos = sv.find( 1UL );

      if( pos != sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Non-existing element could be found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the lowerBound member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the lowerBound member function of SparseSubvector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   test_ = "SparseSubvector::lowerBound()";

   initialize();

   SVT sv = subvector( vec_, 0UL, 3UL );

   // Determining the lower bound for index 0
   {
      SVT::Iterator pos = sv.lowerBound( 0UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 1
   {
      SVT::Iterator pos = sv.lowerBound( 1UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 2
   {
      SVT::Iterator pos = sv.lowerBound( 2UL );

      if( pos != sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the upperBound member function of SparseSubvector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the upperBound member function of SparseSubvector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   test_ = "SparseSubvector::upperBound()";

   initialize();

   SVT sv = subvector( vec_, 0UL, 3UL );

   // Determining the upper bound for index 0
   {
      SVT::Iterator pos = sv.upperBound( 0UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 1
   {
      SVT::Iterator pos = sv.upperBound( 1UL );

      if( pos != sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 2
   {
      SVT::Iterator pos = sv.upperBound( 2UL );

      if( pos != sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDefault function with the SparseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDefault function with the SparseSubvector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   test_ = "isDefault() function";

   initialize();

   // isDefault with default vector
   {
      VT vec( 8UL );
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
/*!\brief Test of the isnan function with the SparseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isnan function with the SparseSubvector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsNan()
{
   test_ = "isnan() function";

   typedef blaze::CompressedVector<float,blaze::columnVector>  VectorType;
   typedef blaze::SparseSubvector<VectorType>                  SubvectorType;

   VectorType vec( 9UL, 4UL );
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
/*!\brief Test of the min function with the SparseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the min function used with the SparseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMinimum()
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
/*!\brief Test of the max function with the SparseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the max function used with the SparseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMaximum()
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
/*!\brief Test of the subvector function with the SparseSubvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subvector function used with the SparseSubvector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubvector()
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

      if( sv2.begin()->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << sv2.begin()->value() << "\n"
             << "   Expected result: -2\n";
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
void ClassTest::initialize()
{
   // Initializing the compressed row vector
   vec_.reset();
   vec_[1] =  1;
   vec_[3] = -2;
   vec_[4] = -3;
   vec_[6] =  4;
}
//*************************************************************************************************

} // namespace sparsesubvector

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
   std::cout << "   Running SparseSubvector class test..." << std::endl;

   try
   {
      RUN_SPARSESUBVECTOR_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SparseSubvector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
