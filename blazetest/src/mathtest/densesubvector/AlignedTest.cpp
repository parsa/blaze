//=================================================================================================
/*!
//  \file src/mathtest/densesubvector/AlignedTest.cpp
//  \brief Source file for the aligned DenseSubvector class test
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
#include <blazetest/mathtest/densesubvector/AlignedTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace mathtest {

namespace densesubvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the aligned DenseSubvector class test.
//
// \exception std::runtime_error Operation error detected.
*/
AlignedTest::AlignedTest()
   : vec1_( 64UL )
   , vec2_( 64UL )
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
void AlignedTest::testConstructors()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "DenseSubvector constructor";

   initialize();

   const size_t alignment = blaze::AlignmentTrait<int>::value;

   for( size_t start=0UL; start<vec1_.size(); start+=alignment ) {
      for( size_t maxsize=0UL; ; maxsize+=alignment )
      {
         const size_t size( blaze::min( maxsize, vec1_.size()-start ) );

         const ASVT sv1 = subvector<aligned>  ( vec1_, start, size );
         const USVT sv2 = subvector<unaligned>( vec2_, start, size );

         if( sv1 != sv2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of dense subvector failed\n"
                << " Details:\n"
                << "   Start = " << start << "\n"
                << "   Size  = " << size  << "\n"
                << "   Subvector:\n" << sv1 << "\n"
                << "   Reference:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( start+maxsize > vec1_.size() ) break;
      }
   }

   try {
      ASVT sv = subvector<aligned>( vec1_, 8UL, 64UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      ASVT sv = subvector<aligned>( vec1_, 80UL, 0UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      ASVT sv = subvector<aligned>( vec1_, 7UL, 16UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of unaligned subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      ASVT sv = subvector<aligned>( vec1_, 8UL, 13UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of unaligned subvector succeeded\n"
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
void AlignedTest::testAssignment()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   //=====================================================================================
   // Homogeneous assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector homogeneous assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      sv1 = 12;
      sv2 = 12;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector copy assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2, 8UL, 16UL );
      sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "DenseSubvector copy assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      sv1 = subvector( vec1_, 24UL, 16UL );
      sv2 = subvector( vec2_, 24UL, 16UL );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector dense vector assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 16UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 = vec;
      sv2 = vec;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector sparse vector assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 16UL );
      randomize( vec, 6UL, int(randmin), int(randmax) );

      sv1 = vec;
      sv2 = vec;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testAddAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   //=====================================================================================
   // DenseSubvector addition assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector addition assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2, 8UL, 16UL );
      sv1 += subvector<aligned>  ( vec1_, 8UL, 16UL );
      sv2 += subvector<unaligned>( vec2_, 8UL, 16UL );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "DenseSubvector addition assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      sv1 += subvector<aligned>  ( vec1_, 24UL, 16UL );
      sv2 += subvector<unaligned>( vec2_, 24UL, 16UL );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector dense vector addition assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 16UL, 0 );
      randomize( vec, int(randmin), int(randmax) );

      sv1 += vec;
      sv2 += vec;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector sparse vector addition assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 16UL );
      randomize( vec, 6UL, int(randmin), int(randmax) );

      sv1 += vec;
      sv2 += vec;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testSubAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   //=====================================================================================
   // DenseSubvector subtraction assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector subtraction assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2, 8UL, 16UL );
      sv1 -= subvector<aligned>  ( vec1_, 24UL, 16UL );
      sv2 -= subvector<unaligned>( vec2_, 24UL, 16UL );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "DenseSubvector subtraction assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      sv1 -= subvector<aligned>  ( vec1_, 24UL, 16UL );
      sv2 -= subvector<unaligned>( vec2_, 24UL, 16UL );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector dense vector subtraction assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 16UL, 0 );
      randomize( vec, int(randmin), int(randmax) );

      sv1 -= vec;
      sv2 -= vec;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector sparse vector subtraction assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 16UL );
      randomize( vec, 6UL, int(randmin), int(randmax) );

      sv1 -= vec;
      sv2 -= vec;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testMultAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   //=====================================================================================
   // DenseSubvector multiplication assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector multiplication assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2, 8UL, 16UL );
      sv1 *= subvector<aligned>  ( vec1_, 24UL, 16UL );
      sv2 *= subvector<unaligned>( vec2_, 24UL, 16UL );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "DenseSubvector multiplication assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      sv1 *= subvector<aligned>  ( vec1_, 24UL, 16UL );
      sv2 *= subvector<unaligned>( vec2_, 24UL, 16UL );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector dense vector multiplication assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 16UL, 0 );
      randomize( vec, int(randmin), int(randmax) );

      sv1 *= vec;
      sv2 *= vec;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector sparse vector multiplication assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 16UL );
      randomize( vec, 6UL, int(randmin), int(randmax) );

      sv1 -= vec;
      sv2 -= vec;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector scalar multiplication assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      sv1 *= 3;
      sv2 *= 3;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testDivAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   //=====================================================================================
   // Scalar division assignment
   //=====================================================================================

   {
      test_ = "DenseSubvector scalar division assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      sv1 /= 0.5;
      sv2 /= 0.5;

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testSubscript()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "DenseSubvector::operator[]";

   initialize();

   ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
   USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

   // Writing the first element
   sv1[1] = 9;
   sv2[1] = 9;

   checkSize( sv1, 16UL );
   checkSize( sv2, 16UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Writing the second element
   sv1[2] = 0;
   sv2[2] = 0;

   checkSize( sv1, 16UL );
   checkSize( sv2, 16UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Writing the third element
   sv1[3] = -8;
   sv2[3] = -8;

   checkSize( sv1, 16UL );
   checkSize( sv2, 16UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testIterator()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   initialize();

   // Counting the number of elements in first half of the vector
   {
      test_ = "Iterator subtraction";

      ASVT sv = subvector<aligned>( vec1_, 0UL, 16UL );
      const size_t number( sv.end() - sv.begin() );

      if( number != 16UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 16\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in second half of the vector
   {
      test_ = "Iterator subtraction";

      ASVT sv = subvector<aligned>( vec1_, 16UL, 48UL );
      const size_t number( sv.end() - sv.begin() );

      if( number != 48UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 48\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      ASVT sv = subvector<aligned>( vec1_, 8UL, 8UL );
      ASVT::ConstIterator it ( sv.cbegin() );
      ASVT::ConstIterator end( sv.cend() );

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || *it != sv[1] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      --it;

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it == end || *it != sv[1] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it--;

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it += 2UL;

      if( it == end || *it != sv[2] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator addition assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it -= 2UL;

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator subtraction assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it + 3UL;

      if( it == end || *it != sv[3] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar addition failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it - 3UL;

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar subtraction failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = 8UL + it;

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

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      int value = 6;

      ASVT::Iterator it1( sv1.begin() );
      USVT::Iterator it2( sv2.begin() );

      for( ; it1!=sv1.end(); ++it1, ++it2 ) {
         *it1 = value;
         *it2 = value;
         ++value;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      int value = 6;

      ASVT::Iterator it1( sv1.begin() );
      USVT::Iterator it2( sv2.begin() );

      for( ; it1!=sv1.end(); ++it1, ++it2 ) {
         *it1 += value;
         *it2 += value;
         ++value;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      int value = 6;

      ASVT::Iterator it1( sv1.begin() );
      USVT::Iterator it2( sv2.begin() );

      for( ; it1!=sv1.end(); ++it1, ++it2 ) {
         *it1 -= value;
         *it2 -= value;
         ++value;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );
      int value = 1;

      ASVT::Iterator it1( sv1.begin() );
      USVT::Iterator it2( sv2.begin() );

      for( ; it1!=sv1.end(); ++it1, ++it2 ) {
         *it1 *= value;
         *it2 *= value;
         ++value;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

      ASVT::Iterator it1( sv1.begin() );
      USVT::Iterator it2( sv2.begin() );

      for( ; it1!=sv1.end(); ++it1, ++it2 ) {
         *it1 /= 2;
         *it2 /= 2;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testNonZeros()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "DenseSubvector::nonZeros()";

   initialize();

   // Initialization check
   ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
   USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

   checkSize( sv1, 16UL );
   checkSize( sv2, 16UL );

   if( sv1.nonZeros() != sv2.nonZeros() ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense subvector
   sv1[3] = 0;
   sv2[3] = 0;

   checkSize( sv1, 16UL );
   checkSize( sv2, 16UL );

   if( sv1.nonZeros() != sv2.nonZeros() ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense vector
   vec1_[9UL] = 5;
   vec2_[9UL] = 5;

   checkSize( sv1, 16UL );
   checkSize( sv2, 16UL );

   if( sv1.nonZeros() != sv2.nonZeros() ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testReset()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "DenseSubvector::reset()";

   initialize();

   // Resetting the range [0,15]
   {
      ASVT sv1 = subvector<aligned>  ( vec1_, 0UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 0UL, 16UL );
      sv1.reset();
      sv2.reset();

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [0,15] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting the range [16,63]
   {
      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 48UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 48UL );
      sv1.reset();
      sv2.reset();

      checkSize( sv1, 48UL );
      checkSize( sv2, 48UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [16,63] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testScale()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "DenseSubvector::scale()";

   initialize();

   ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 16UL );
   USVT sv2 = subvector<unaligned>( vec2_, 8UL, 16UL );

   // Integral scaling of the subvector in the range [8,23]
   sv1.scale( 3 );
   sv2.scale( 3 );

   checkSize( sv1, 16UL );
   checkSize( sv2, 16UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Integral scale operation of range [8,23] failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Floating point scaling of the subvector in the range [8,23]
   sv1.scale( 0.5 );
   sv2.scale( 0.5 );

   checkSize( sv1, 16UL );
   checkSize( sv2, 16UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Floating point scale operation of range [8,23] failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
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
void AlignedTest::testIsDefault()
{
   using blaze::subvector;
   using blaze::aligned;


   test_ = "isDefault() function";

   initialize();

   // isDefault with default vector
   {
      VT vec( 64UL, 0 );
      ASVT sv = subvector<aligned>( vec, 8UL, 16UL );

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
      ASVT sv = subvector<aligned>( vec1_, 8UL, 16UL );

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
void AlignedTest::testIsNan()
{
   using blaze::subvector;
   using blaze::aligned;


   test_ = "isnan() function";

   typedef blaze::DynamicVector<float,blaze::rowVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType,aligned>     SubvectorType;

   VectorType vec( vec1_ );
   subvector<aligned>( vec, 0UL, 32UL ) = 0;

   // isnan with empty 32-dimensional subvector
   {
      SubvectorType sv = subvector<aligned>( vec, 0UL, 32UL );

      checkSize    ( sv, 32UL );
      checkNonZeros( sv,  0UL );

      if( blaze::isnan( sv ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with fully filled 32-dimensional subvector
   {
      SubvectorType sv = subvector<aligned>( vec, 32UL, 32UL );

      checkSize( sv, 32UL );

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
void AlignedTest::testMinimum()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "min() function";

   initialize();

   // Computing the minimum of the in the range [0,15]
   {
      const int minimum1 = min( subvector<aligned>  ( vec1_, 0UL, 16UL ) );
      const int minimum2 = min( subvector<unaligned>( vec2_, 0UL, 16UL ) );

      if( minimum1 != minimum2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Minimum computation for range [0,15] failed\n"
             << " Details:\n"
             << "   Result: " << minimum1 << "\n"
             << "   Expected result: " << minimum2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the minimum of the in the range [16,31]
   {
      const int minimum1 = min( subvector<aligned>  ( vec1_, 16UL, 16UL ) );
      const int minimum2 = min( subvector<unaligned>( vec2_, 16UL, 16UL ) );

      if( minimum1 != minimum2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Minimum computation for range [16,31] failed\n"
             << " Details:\n"
             << "   Result: " << minimum1 << "\n"
             << "   Expected result: " << minimum2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the minimum of the in the range [32,47]
   {
      const int minimum1 = min( subvector<aligned>  ( vec1_, 32UL, 16UL ) );
      const int minimum2 = min( subvector<unaligned>( vec2_, 32UL, 16UL ) );

      if( minimum1 != minimum2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Minimum computation for range [32,47] failed\n"
             << " Details:\n"
             << "   Result: " << minimum1 << "\n"
             << "   Expected result: " << minimum2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the minimum of the in the range [48,63]
   {
      const int minimum1 = min( subvector<aligned>  ( vec1_, 48UL, 16UL ) );
      const int minimum2 = min( subvector<unaligned>( vec2_, 48UL, 16UL ) );

      if( minimum1 != minimum2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Minimum computation for range [48,63] failed\n"
             << " Details:\n"
             << "   Result: " << minimum1 << "\n"
             << "   Expected result: " << minimum2 << "\n";
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
void AlignedTest::testMaximum()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "max() function";

   initialize();

   // Computing the maximum of the in the range [0,15]
   {
      const int maximum1 = max( subvector<aligned>  ( vec1_, 0UL, 16UL ) );
      const int maximum2 = max( subvector<unaligned>( vec2_, 0UL, 16UL ) );

      if( maximum1 != maximum2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Maximum computation for range [0,15] failed\n"
             << " Details:\n"
             << "   Result: " << maximum1 << "\n"
             << "   Expected result: " << maximum2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the maximum of the in the range [16,31]
   {
      const int maximum1 = max( subvector<aligned>  ( vec1_, 16UL, 16UL ) );
      const int maximum2 = max( subvector<unaligned>( vec2_, 16UL, 16UL ) );

      if( maximum1 != maximum2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Maximum computation for range [16,31] failed\n"
             << " Details:\n"
             << "   Result: " << maximum1 << "\n"
             << "   Expected result: " << maximum2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the maximum of the in the range [32,47]
   {
      const int maximum1 = max( subvector<aligned>  ( vec1_, 32UL, 16UL ) );
      const int maximum2 = max( subvector<unaligned>( vec2_, 32UL, 16UL ) );

      if( maximum1 != maximum2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Maximum computation for range [32,47] failed\n"
             << " Details:\n"
             << "   Result: " << maximum1 << "\n"
             << "   Expected result: " << maximum2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Computing the maximum of the in the range [48,63]
   {
      const int maximum1 = max( subvector<aligned>  ( vec1_, 48UL, 16UL ) );
      const int maximum2 = max( subvector<unaligned>( vec2_, 48UL, 16UL ) );

      if( maximum1 != maximum2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Maximum computation for range [48,63] failed\n"
             << " Details:\n"
             << "   Result: " << maximum1 << "\n"
             << "   Expected result: " << maximum2 << "\n";
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
void AlignedTest::testSubvector()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "subvector() function";

   initialize();

   {
      ASVT sv1 = subvector<aligned>  ( vec1_, 8UL, 32UL );
      ASVT sv2 = subvector<aligned>  ( sv1  , 8UL, 16UL );
      USVT sv3 = subvector<unaligned>( vec2_, 8UL, 32UL );
      USVT sv4 = subvector<unaligned>( sv3  , 8UL, 16UL );

      if( sv2 != sv4 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subvector function failed\n"
             << " Details:\n"
             << "   Result:\n" << sv2 << "\n"
             << "   Expected result:\n" << sv4 << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv2[1] != sv4[1] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << sv2[1] << "\n"
             << "   Expected result: " << sv4[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( *sv2.begin() != *sv4.begin() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *sv2.begin() << "\n"
             << "   Expected result: " << *sv4.begin() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ASVT sv1 = subvector<aligned>( vec1_,  8UL, 32UL );
      ASVT sv2 = subvector<aligned>( sv1  , 32UL,  8UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      ASVT sv1 = subvector<aligned>( vec1_, 8UL, 32UL );
      ASVT sv2 = subvector<aligned>( sv1  , 8UL, 32UL );

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
void AlignedTest::initialize()
{
   // Initializing the dynamic row vectors
   randomize( vec1_, int(randmin), int(randmax) );
   vec2_ = vec1_;
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
   std::cout << "   Running aligned DenseSubvector class test..." << std::endl;

   try
   {
      RUN_DENSESUBVECTOR_ALIGNED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aligned DenseSubvector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
