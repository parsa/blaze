//=================================================================================================
/*!
//  \file src/mathtest/compressedvector/ClassTest.cpp
//  \brief Source file for the CompressedVector class test
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
#include <blaze/math/shims/Equal.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/compressedvector/ClassTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace mathtest {

namespace compressedvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the CompressedVector class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testSubscript();
   testNonZeros();
   testReset();
   testClear();
   testAppend();
   testInsert();
   testErase();
   testResize();
   testReserve();
   testScale();
   testSwap();
   testFind();
   testLowerBound();
   testUpperBound();
   testIsDefault();
   testIsNan();
   testLength();
   testNormalize();
   testMinimum();
   testMaximum();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the CompressedVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the CompressedVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "CompressedVector default constructor";

      blaze::CompressedVector<int,blaze::rowVector> vec;

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // Size constructor
   //=====================================================================================

   {
      test_ = "CompressedVector size constructor (size 0)";

      blaze::CompressedVector<int,blaze::rowVector> vec( 0UL );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "CompressedVector size constructor (size 5)";

      blaze::CompressedVector<int,blaze::rowVector> vec( 5UL );

      checkSize    ( vec, 5UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // Size/Non-zeros constructor
   //=====================================================================================

   {
      test_ = "CompressedVector size/non-zeros constructor (size 0)";

      blaze::CompressedVector<int,blaze::rowVector> vec( 0UL, 3UL );

      checkSize    ( vec, 0UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "CompressedVector size/non-zeros constructor (size 7)";

      blaze::CompressedVector<int,blaze::rowVector> vec( 7UL, 3UL );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "CompressedVector copy constructor (size 0)";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 0UL, 3UL );
      blaze::CompressedVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 0UL );
      checkNonZeros( vec2, 0UL );
   }

   {
      test_ = "CompressedVector copy constructor (size 7)";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 7UL, 3UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[3] = 4;
      blaze::CompressedVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 7UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the CompressedVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "CompressedVector copy assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 7UL, 3UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[3] = 4;
      blaze::CompressedVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 7UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CompressedVector copy assignment stress test";

      typedef blaze::CompressedVector<int,blaze::rowVector>  RandomVectorType;

      blaze::CompressedVector<int,blaze::rowVector> vec1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 20UL ) );
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( size, min, max ) );

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
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1( 5UL, 0 );
      vec1[0] = 10;
      vec1[1] = 11;
      vec1[2] = 12;
      vec1[4] = 13;
      blaze::CompressedVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 10 || vec2[1] != 11 || vec2[2] != 12 || vec2[3] != 0 || vec2[4] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 10 11 12 0 13 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CompressedVector dense vector assignment stress test";

      typedef blaze::DynamicVector<int,blaze::rowVector>  RandomVectorType;

      blaze::CompressedVector<int,blaze::rowVector> vec1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 20UL ) );
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( size, min, max ) );

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
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1( 7UL, 3UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[3] = 4;
      blaze::CompressedVector<int,blaze::rowVector> vec2;
      vec2 = trans( vec1 );

      checkSize    ( vec2, 7UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CompressedVector sparse vector assignment stress test";

      typedef blaze::CompressedVector<unsigned int,blaze::rowVector>  RandomVectorType;

      blaze::CompressedVector<int,blaze::rowVector> vec1;
      const unsigned int min( randmin );
      const unsigned int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 20UL ) );
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( size, min, max ) );

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
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector addition assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1( 5UL, 0 );
      vec1[0] = 10;
      vec1[1] = 11;
      vec1[2] = 12;
      vec1[4] = 13;
      blaze::CompressedVector<int,blaze::rowVector> vec2( 5UL, 3UL );
      vec2[0] = 1;
      vec2[1] = 2;
      vec2[3] = 4;

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 11 || vec2[1] != 13 || vec2[2] != 12 || vec2[3] != 4 || vec2[4] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 11 13 12 4 13 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector addition assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1( 5UL, 3UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[3] = 4;
      blaze::CompressedVector<int,blaze::rowVector> vec2( 5UL, 2UL );
      vec2[1] = 5;
      vec2[2] = 6;

      vec2 += trans( vec1 );

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 7 || vec2[2] != 6 || vec2[3] != 4 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 7 6 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector subtraction assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1( 5UL, 0 );
      vec1[0] = 10;
      vec1[1] = 11;
      vec1[2] = 12;
      vec1[4] = 13;
      blaze::CompressedVector<int,blaze::rowVector> vec2( 5UL, 3UL );
      vec2[0] = 1;
      vec2[1] = 2;
      vec2[3] = 4;

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != -9 || vec2[1] != -9 || vec2[2] != -12 || vec2[3] != 4 || vec2[4] != -13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( -9 -9 -12 4 -13 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector subtraction assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1( 5UL, 3UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[3] = 4;
      blaze::CompressedVector<int,blaze::rowVector> vec2( 5UL, 2UL );
      vec2[1] = 5;
      vec2[2] = 6;

      vec2 -= trans( vec1 );

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != -1 || vec2[1] != 3 || vec2[2] != 6 || vec2[3] != -4 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( -1 3 6 -4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector multiplication assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1( 5UL, 0 );
      vec1[0] = 10;
      vec1[1] = 11;
      vec1[2] = 12;
      vec1[4] = 13;
      blaze::CompressedVector<int,blaze::rowVector> vec2( 5UL, 3UL );
      vec2[0] = 1;
      vec2[1] = 2;
      vec2[3] = 4;

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 10 || vec2[1] != 22 || vec2[2] != 0 || vec2[3] != 0 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 10 22 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector multiplication assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1( 5UL, 3UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[3] = 4;
      blaze::CompressedVector<int,blaze::rowVector> vec2( 5UL, 2UL );
      vec2[1] = 5;
      vec2[2] = 6;

      vec2 *= trans( vec1 );

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 1UL );

      if( vec2[0] != 0 || vec2[1] != 10 || vec2[2] != 0 || vec2[3] != 0 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 10 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "CompressedVector scalar multiplication assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1( 5UL, 3UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[3] = 4;

      vec1 *= 2;

      checkSize    ( vec1, 5UL );
      checkNonZeros( vec1, 3UL );

      if( vec1[0] != 2 || vec1[1] != 4 || vec1[2] != 0 || vec1[3] != 8 || vec1[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 2 4 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testDivAssign()
{
   //=====================================================================================
   // Scalar division assignment
   //=====================================================================================

   {
      test_ = "CompressedVector scalar division assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1( 5UL, 3UL );
      vec1[0] = 2;
      vec1[1] = 4;
      vec1[3] = 8;

      vec1 /= 2;

      checkSize    ( vec1, 5UL );
      checkNonZeros( vec1, 3UL );

      if( vec1[0] != 1 || vec1[1] != 2 || vec1[2] != 0 || vec1[3] != 4 || vec1[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the CompressedVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSubscript()
{
   test_ = "CompressedVector::operator[]";

   // Adding the first element
   blaze::CompressedVector<int,blaze::rowVector> vec( 7UL, 3UL );
   vec[2] = 1;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 1UL );

   if( vec[2] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Adding the second element
   vec[5] = 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 2UL );

   if( vec[2] != 1 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Adding the third element
   vec[3] = 3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   if( vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Adding the fourth element
   vec[0] = 4;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the nonZeros member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   test_ = "CompressedVector::nonZeros()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 7UL, 3UL );

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 0UL );

   // Adding two non-zero elements
   vec[2] = 1;
   vec[5] = 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 2UL );

   // Adding a third non-zero elements
   vec[3] = 0;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 2UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reset member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   test_ = "CompressedVector::reset()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 11UL, 4UL );
   vec[1] = 1;
   vec[3] = 2;
   vec[7] = 3;
   vec[9] = 4;

   checkSize    ( vec, 11UL );
   checkCapacity( vec,  4UL );
   checkNonZeros( vec,  4UL );

   if( vec[1] != 1 || vec[3] != 2 || vec[7] != 3 || vec[9] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 1 0 3 0 0 0 3 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resetting the vector
   vec.reset();

   checkSize    ( vec, 11UL );
   checkNonZeros( vec,  0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the clear member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the clear member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   test_ = "CompressedVector::clear()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 9UL, 3UL );
   vec[0] = 1;
   vec[7] = 2;
   vec[8] = 3;

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   if( vec[0] != 1 || vec[7] != 2 || vec[8] != 3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 0 0 0 0 0 0 2 3 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Clearing the vector
   vec.clear();

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the append member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the append member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAppend()
{
   test_ = "CompressedVector::append()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 9UL, 4UL );

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 0UL );

   // Appending one non-zero element
   vec.append( 1UL, 1 );

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 1UL );

   if( vec[1] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 1 0 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Appending three more non-zero elements
   vec.append( 3UL, 2 );
   vec.append( 4UL, 3 );
   vec.append( 8UL, 4 );

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec[1] != 1 || vec[3] != 2 || vec[4] != 3 || vec[8] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Append operation failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 1 0 2 3 0 0 0 4 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the insert member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the insert member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testInsert()
{
   test_ = "CompressedVector::insert()";

   typedef blaze::CompressedVector<int,blaze::rowVector>::Iterator  Iterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 7UL );

   checkSize    ( vec, 7UL );
   checkNonZeros( vec, 0UL );

   // Inserting a non-zero element
   {
      Iterator pos = vec.insert( 4UL, 1 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      if( pos->value() != 1 || pos->index() != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 1\n"
             << "   Expected index: 4\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[4] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a second non-zero element
   {
      Iterator pos = vec.insert( 6UL, 2 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( pos->value() != 2 || pos->index() != 6UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 2\n"
             << "   Expected index: 6\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a third non-zero element
   {
      Iterator pos = vec.insert( 2UL, 3 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( pos->value() != 3 || pos->index() != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 3\n"
             << "   Expected index: 2\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[2] != 3 || vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 3 0 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a fourth non-zero element
   {
      Iterator pos = vec.insert( 3UL, 4 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( pos->value() != 4 || pos->index() != 3UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 4\n"
             << "   Expected index: 3\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[2] != 3 || vec[3] != 4 || vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 3 4 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Trying to insert an already existing element
   try {
      vec.insert( 3UL, 5 );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Inserting an existing element succeeded\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 3 4 1 0 2 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the erase member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the erase member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testErase()
{
   //=====================================================================================
   // Index-based erase function
   //=====================================================================================

   {
      test_ = "CompressedVector::erase( size_t )";

      // Initialization check
      blaze::CompressedVector<int,blaze::rowVector> vec( 9UL, 5UL );
      vec[0] = 1;
      vec[2] = 2;
      vec[5] = 3;
      vec[7] = 4;
      vec[8] = 5;

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 1 || vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at index 0
      vec.erase( size_t(0) );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 4UL );

      if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at index 8
      vec.erase( 8UL );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at index 5
      vec.erase( 5UL );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 2UL );

      if( vec[2] != 2 || vec[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      vec.erase( 1UL );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 2UL );

      if( vec[2] != 2 || vec[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Iterator-based erase() function
   //=====================================================================================

   {
      test_ = "CompressedVector::erase( Iterator )";

      typedef blaze::CompressedVector<int,blaze::rowVector>  VectorType;
      typedef VectorType::Iterator  Iterator;

      // Initialization check
      VectorType vec( 9UL, 5UL );
      vec[0] = 1;
      vec[2] = 2;
      vec[5] = 3;
      vec[7] = 4;
      vec[8] = 5;

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 1 || vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at index 0
      {
         Iterator pos = vec.erase( vec.find( 0UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 4UL );

         if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 3 0 4 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 2 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at index 8
      {
         Iterator pos = vec.erase( vec.find( 8UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 3UL );

         if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at index 5
      {
         Iterator pos = vec.erase( vec.find( 5UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 2UL );

         if( vec[2] != 2 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 4 || pos->index() != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 7\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = vec.erase( vec.find( 1UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 2UL );

         if( vec[2] != 2 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Iterator-range-based erase() function
   //=====================================================================================

   {
      test_ = "CompressedVector::erase( Iterator, Iterator )";

      typedef blaze::CompressedVector<int,blaze::rowVector>  VectorType;
      typedef VectorType::Iterator  Iterator;

      // Initialization check
      VectorType vec( 9UL, 5UL );
      vec[0] = 1;
      vec[2] = 2;
      vec[5] = 3;
      vec[7] = 4;
      vec[8] = 5;

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 1 || vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the range from index 0 to index 2
      {
         Iterator pos = vec.erase( vec.find( 0UL ), vec.find( 2UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 4UL );

         if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 3 0 4 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 2 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the range from index 8 to the end
      {
         Iterator pos = vec.erase( vec.find( 8UL ), vec.end() );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 3UL );

         if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the range from index 5 to index 7
      {
         Iterator pos = vec.erase( vec.find( 5UL ), vec.find( 7UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 2UL );

         if( vec[2] != 2 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 4 || pos->index() != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 7\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = vec.erase( vec.find( 2UL ), vec.find( 2UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 2UL );

         if( vec[2] != 2 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != vec.find( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the resize member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the resize member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   test_ = "CompressedVector::resize()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 0
   vec.resize( 0UL );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 5
   vec.resize( 5UL );

   checkSize    ( vec, 5UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 9 and preserving the elements
   vec[0] = 1;
   vec[2] = 2;
   vec[4] = 3;
   vec.resize( 9UL, true );

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   if( vec[0] != 1 || vec[2] != 2 || vec[4] != 3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 0 2 0 3 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 2 and preserving the elements
   vec.resize( 2UL, true );

   checkSize    ( vec, 2UL );
   checkCapacity( vec, 1UL );
   checkNonZeros( vec, 1UL );

   if( vec[0] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 1
   vec.resize( 1UL );

   checkSize( vec, 1UL );

   // Resizing to 0
   vec.resize( 0UL );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reserve member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reserve member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReserve()
{
   test_ = "CompressedVector::reserve()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Increasing the capacity of the vector
   vec.reserve( 10UL );

   checkSize    ( vec,  0UL );
   checkCapacity( vec, 10UL );
   checkNonZeros( vec,  0UL );

   // Further increasing the capacity of the vector
   vec.reserve( 20UL );

   checkSize    ( vec,  0UL );
   checkCapacity( vec, 20UL );
   checkNonZeros( vec,  0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of CompressedVector.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScale()
{
   test_ = "CompressedVector::scale()";

   {
      // Initialization check
      blaze::CompressedVector<int,blaze::rowVector> vec( 6UL, 3UL );
      vec[1] = 1;
      vec[3] = 2;
      vec[5] = 3;

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[1] != 1 || vec[3] != 2 || vec[5] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 1 0 2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the vector
      vec.scale( 2 );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[1] != 2 || vec[3] != 4 || vec[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 2 0 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the vector
      vec.scale( 0.5 );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[1] != 1 || vec[3] != 2 || vec[5] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 1 0 2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      using blaze::complex;

      blaze::CompressedVector<complex<float>,blaze::rowVector> vec( 2UL, 2UL );
      vec[0] = complex<float>( 1.0F, 0.0F );
      vec[1] = complex<float>( 2.0F, 0.0F );
      vec.scale( complex<float>( 3.0F, 0.0F ) );

      checkSize    ( vec, 2UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != complex<float>( 3.0F, 0.0F ) || vec[1] != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( (3,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the swap functionality of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the swap function of the CompressedVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   test_ = "CompressedVector swap";

   blaze::CompressedVector<int,blaze::rowVector> vec1( 12UL, 4UL );
   vec1[ 1] = 1;
   vec1[ 4] = 2;
   vec1[ 7] = 3;
   vec1[10] = 4;

   blaze::CompressedVector<int,blaze::rowVector> vec2( 5UL, 2UL );
   vec2[0] = 4;
   vec2[4] = 2;

   swap( vec1, vec2 );

   checkSize    ( vec1, 5UL );
   checkCapacity( vec1, 2UL );
   checkNonZeros( vec1, 2UL );

   if( vec1[0] != 4 || vec1[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 4 0 0 0 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   checkSize    ( vec2, 12UL );
   checkCapacity( vec2,  4UL );
   checkNonZeros( vec2,  4UL );

   if( vec2[1] != 1 || vec2[4] != 2 || vec2[7] != 3 || vec2[10] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 0 1 0 0 2 0 0 3 0 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the find member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the find member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   test_ = "CompressedVector::find()";

   typedef blaze::CompressedVector<int,blaze::rowVector>::ConstIterator  ConstIterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 3UL );
   vec[0] = 1;
   vec[2] = 2;
   vec[7] = 3;

   checkSize    ( vec, 8UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   // Searching for the first element
   {
      ConstIterator pos( vec.find( 0UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current vector:\n" << vec << "\n";
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
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the second element
   {
      ConstIterator pos( vec.find( 2UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the third element
   {
      ConstIterator pos( vec.find( 7UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for a non-existing non-zero element
   {
      ConstIterator pos( vec.find( 5UL ) );

      if( pos != vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Non-existing element could be found\n"
             << " Details:\n"
             << "   Required index = 5\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 0\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the lowerBound member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the lowerBound member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   test_ = "CompressedVector::lowerBound()";

   typedef blaze::CompressedVector<int,blaze::rowVector>::ConstIterator  ConstIterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 3UL );
   vec[0] = 1;
   vec[2] = 2;
   vec[7] = 3;

   checkSize    ( vec, 8UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   // Determining the lower bound for index 0
   {
      ConstIterator pos( vec.lowerBound( 0UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current vector:\n" << vec << "\n";
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
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 1
   {
      ConstIterator pos( vec.lowerBound( 1UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 2
   {
      ConstIterator pos( vec.lowerBound( 2UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 3
   {
      ConstIterator pos( vec.lowerBound( 3UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 7
   {
      ConstIterator pos( vec.lowerBound( 7UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the upperBound member function of CompressedVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the upperBound member function of CompressedVector. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   test_ = "CompressedVector::upperBound()";

   typedef blaze::CompressedVector<int,blaze::rowVector>::ConstIterator  ConstIterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 3UL );
   vec[0] = 1;
   vec[2] = 2;
   vec[7] = 3;

   checkSize    ( vec, 8UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   // Determining the upper bound for index 0
   {
      ConstIterator pos( vec.upperBound( 0UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 1
   {
      ConstIterator pos( vec.upperBound( 1UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 2
   {
      ConstIterator pos( vec.upperBound( 2UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 3
   {
      ConstIterator pos( vec.upperBound( 3UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 7
   {
      ConstIterator pos( vec.upperBound( 7UL ) );

      if( pos != vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDefault function with the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDefault function with the CompressedVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   test_ = "isDefault() function";

   // isDefault with vector of size 0
   {
      blaze::CompressedVector<int,blaze::rowVector> vec;

      if( isDefault( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with default vector
   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL );

      if( isDefault( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with non-default vector
   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 1UL );
      vec[1] = 1;

      if( isDefault( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 2 function with the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isnan function with the CompressedVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsNan()
{
   test_ = "isnan() function";

   // isnan with 0-dimensional vector
   {
      blaze::CompressedVector<float,blaze::rowVector> vec;

      if( blaze::isnan( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with empty 9-dimensional vector
   {
      blaze::CompressedVector<float,blaze::rowVector> vec( 9UL );

      if( blaze::isnan( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with filled 9-dimensional vector
   {
      blaze::CompressedVector<float,blaze::rowVector> vec( 9UL );
      vec[3] =  1.0F;
      vec[4] = -2.0F;
      vec[6] =  3.0F;
      vec[8] =  4.0F;

      if( blaze::isnan( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the length and sqrLength functions with the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the length and sqrLength functions with the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLength()
{
   test_ = "length() and sqrLength() functions";

   {
      // Initialization check
      blaze::CompressedVector<double,blaze::rowVector> vec;

      // Computing the vector length
      const double len( length( vec ) );

      if( !blaze::equal( len, 0.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Length computation failed\n"
             << " Details:\n"
             << "   Result: " << len << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }

      // Computing the vector square length
      const double sqrlen( sqrLength( vec ) );

      if( !blaze::equal( sqrLength( vec ), 0.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Square length computation failed\n"
             << " Details:\n"
             << "   Result: " << sqrlen << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      // Initialization check
      blaze::CompressedVector<double,blaze::rowVector> vec( 5UL );

      // Computing the vector length
      const double len( length( vec ) );

      if( !blaze::equal( len, 0.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Length computation failed\n"
             << " Details:\n"
             << "   Result: " << len << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }

      // Computing the vector square length
      const double sqrlen( sqrLength( vec ) );

      if( !blaze::equal( sqrLength( vec ), 0.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Square length computation failed\n"
             << " Details:\n"
             << "   Result: " << sqrlen << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      // Initialization check
      blaze::CompressedVector<double,blaze::rowVector> vec( 5UL, 2UL );
      vec[1] = 3.0;
      vec[4] = 4.0;

      // Computing the vector length
      const double len( length( vec ) );

      if( !blaze::equal( len, 5.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Length computation failed\n"
             << " Details:\n"
             << "   Result: " << len << "\n"
             << "   Expected result: 5\n";
         throw std::runtime_error( oss.str() );
      }

      // Computing the vector square length
      const double sqrlen( sqrLength( vec ) );

      if( !blaze::equal( sqrLength( vec ), 25.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Square length computation failed\n"
             << " Details:\n"
             << "   Result: " << sqrlen << "\n"
             << "   Expected result: 25\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the normalize function with the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the normalize function with the CompressedVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNormalize()
{
   test_ = "normalize() function";

   // Initialization check
   blaze::CompressedVector<double,blaze::rowVector> vec( 10UL, 4UL );
   vec[0] = 1.0;
   vec[1] = 2.0;
   vec[2] = 3.0;
   vec[3] = 4.0;

   checkSize    ( vec, 10UL );
   checkCapacity( vec,  4UL );
   checkNonZeros( vec,  4UL );

   if( vec[0] != 1.0 || vec[1] != 2.0 || vec[2] != 3.0 || vec[3] != 4.0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 3 4 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Acquiring normalized vector
   const blaze::CompressedVector<double,blaze::rowVector> normalized( normalize( vec ) );

   if( !blaze::equal( length( normalized ), 1.0 ) ) {
      std::ostringstream oss;
      oss << " Test: CompressedVector::getNormalized()\n"
          << " Error: Normalization failed\n"
          << " Details:\n"
          << "   Result: " << length( normalized ) << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }

   // Normalizing the vector
   vec = normalize( vec );

   if( !blaze::equal( length( vec ), 1.0 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Normalization failed\n"
          << " Details:\n"
          << "   Result: " << length( vec ) << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the min function with the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the min function with the CompressedVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMinimum()
{
   test_ = "min() function";

   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 3UL );
      vec[1] =  1;
      vec[3] =  4;
      vec[7] =  3;

      const int minimum = min( vec );

      if( minimum != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: First computation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 4UL );
      vec[1] =  -4;
      vec[3] =  -2;
      vec[5] =   8;
      vec[7] =  -3;

      const int minimum = min( vec );

      if( minimum != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Second computation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 2UL );
      vec[5] =   8;
      vec[6] =  -3;

      const int minimum = min( vec );

      if( minimum != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Third computation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the max function with the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the max function with the CompressedVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMaximum()
{
   test_ = "max() function";

   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 3UL );
      vec[1] = -1;
      vec[3] = -4;
      vec[7] = -3;

      const int maximum = max( vec );

      if( maximum != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: First computation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 4UL );
      vec[1] =  4;
      vec[3] =  2;
      vec[5] = -8;
      vec[7] =  3;

      const int maximum = max( vec );

      if( maximum != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Second computation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 8UL, 2UL );
      vec[5] = -8;
      vec[6] =  3;

      const int maximum = max( vec );

      if( maximum != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Third computation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 3\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace compressedvector

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
   std::cout << "   Running CompressedVector class test..." << std::endl;

   try
   {
      RUN_COMPRESSEDVECTOR_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during CompressedVector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
