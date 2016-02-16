//=================================================================================================
/*!
//  \file src/mathtest/customvector/AlignedUnpaddedTest.cpp
//  \brief Source file for the aligned/unpadded CustomVector class test
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
#include <blaze/math/shims/Equal.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blaze/util/UniqueArray.h>
#include <blazetest/mathtest/customvector/AlignedUnpaddedTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace mathtest {

namespace customvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the aligned/unpadded CustomVector class test.
//
// \exception std::runtime_error Operation error detected.
*/
AlignedUnpaddedTest::AlignedUnpaddedTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testScaling();
   testSubscript();
   testAt();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testSwap();
   testIsDefault();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the CustomVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the CustomVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testConstructors()
{
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "CustomVector default constructor";

      VT vec;

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // Constructor ( Type*, size_t )
   //=====================================================================================

   {
      test_ = "CustomVector constructor ( Type*, size_t )";

      // Constructing a custom vector of size 10
      {
         blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 10UL ) );
         VT vec( array.get(), 10UL );

         checkSize    ( vec, 10UL );
         checkCapacity( vec, 10UL );
      }

      // Trying to construct a custom vector with invalid array of elements
      try {
         VT vec( NULL, 0UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom vector with a NULL pointer succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom vector with invalid alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 5UL ) );
            VT vec( array.get()+1UL, 4UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom vector with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Constructor ( Type*, size_t, Deleter )
   //=====================================================================================

   {
      test_ = "CustomVector constructor ( Type*, size_t, Deleter )";

      // Constructing a custom vector of size 10
      {
         VT vec( blaze::allocate<int>( 10UL ), 10UL, blaze::Deallocate() );

         checkSize    ( vec, 10UL );
         checkCapacity( vec, 10UL );
      }

      // Trying to construct a custom vector with invalid array of elements
      try {
         VT vec( NULL, 0UL, blaze::Deallocate() );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom vector with a NULL pointer succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom vector with invalid alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 5UL ) );
            VT vec( array.get()+1UL, 4UL, blaze::Deallocate() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom vector with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "CustomVector copy constructor (size 0)";

      VT vec1;
      VT vec2( vec1 );

      checkSize    ( vec2, 0UL );
      checkNonZeros( vec2, 0UL );
   }

   {
      test_ = "CustomVector copy constructor (size 5)";

      VT vec1( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      VT vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec1.data() != vec2.data() ||
          vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomVector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the CustomVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testAssignment()
{
   //=====================================================================================
   // Homogeneous assignment
   //=====================================================================================

   {
      test_ = "CustomVector homogeneous assignment";

      VT vec( blaze::allocate<int>( 3UL ), 3UL, blaze::Deallocate() );
      vec = 2;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Array assignment
   //=====================================================================================

   {
      test_ = "CustomVector array assignment";

      const int array[4] = { 1, 2, 3, 4 };
      VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
      vec = array;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

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


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "CustomVector copy assignment";

      VT vec1( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

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


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "CustomVector dense vector assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      typedef blaze::CustomVector<int,aligned,padded,rowVector>  AlignedPadded;
      AlignedPadded vec1( blaze::allocate<int>( 16UL ), 5UL, 16UL, blaze::Deallocate() );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;

      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

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
      test_ = "CustomVector dense vector assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      typedef blaze::CustomVector<int,unaligned,unpadded,rowVector>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[6UL] );
      UnalignedUnpadded vec1( array.get()+1UL, 5UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;

      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

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


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "CustomVector sparse vector assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL );
      vec1[0] = 1;
      vec1[2] = 2;
      vec1[3] = 3;
      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 0 || vec2[2] != 2 || vec2[3] != 3 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 0 2 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomVector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the CustomVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "CustomVector dense vector addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      typedef blaze::CustomVector<int,aligned,padded,rowVector>  AlignedPadded;
      AlignedPadded vec1( blaze::allocate<int>( 16UL ), 5UL, 16UL, blaze::Deallocate() );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      vec1[3] =  3;
      vec1[4] =  0;

      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CustomVector dense vector addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      typedef blaze::CustomVector<int,unaligned,unpadded,rowVector>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[6UL] );
      UnalignedUnpadded vec1( array.get()+1UL, 5UL );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      vec1[3] =  3;
      vec1[4] =  0;

      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "CustomVector sparse vector addition assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL, 3UL );
      vec1[0] =  1;
      vec1[2] = -2;
      vec1[3] =  3;
      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomVector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the CustomVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "CustomVector dense vector subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      typedef blaze::CustomVector<int,aligned,padded,rowVector>  AlignedPadded;
      AlignedPadded vec1( blaze::allocate<int>( 16UL ), 5UL, 16UL, blaze::Deallocate() );
      vec1[0] = -1;
      vec1[1] =  0;
      vec1[2] =  2;
      vec1[3] = -3;
      vec1[4] =  0;

      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CustomVector dense vector subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      typedef blaze::CustomVector<int,unaligned,unpadded,rowVector>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[6UL] );
      UnalignedUnpadded vec1( array.get()+1UL, 5UL );
      vec1[0] = -1;
      vec1[1] =  0;
      vec1[2] =  2;
      vec1[3] = -3;
      vec1[4] =  0;

      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "CustomVector sparse vector subtraction assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL, 3UL );
      vec1[0] = -1;
      vec1[2] =  2;
      vec1[3] = -3;
      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomVector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the CustomVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "CustomVector dense vector multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      typedef blaze::CustomVector<int,aligned,padded,rowVector>  AlignedPadded;
      AlignedPadded vec1( blaze::allocate<int>( 16UL ), 5UL, 16UL, blaze::Deallocate() );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      vec1[3] =  3;
      vec1[4] =  0;

      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 2UL );

      if( vec2[0] != 0 || vec2[1] != 0 || vec2[2] != -4 || vec2[3] != -18 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 0 -4 -18 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CustomVector dense vector multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      typedef blaze::CustomVector<int,unaligned,unpadded,rowVector>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[6UL] );
      UnalignedUnpadded vec1( array.get()+1UL, 5UL );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      vec1[3] =  3;
      vec1[4] =  0;

      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 2UL );

      if( vec2[0] != 0 || vec2[1] != 0 || vec2[2] != -4 || vec2[3] != -18 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 0 -4 -18 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "CustomVector sparse vector multiplication assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL, 3UL );
      vec1[0] =  1;
      vec1[2] = -2;
      vec1[3] =  3;
      VT vec2( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec2[0] =  0;
      vec2[1] =  4;
      vec2[2] =  2;
      vec2[3] = -6;
      vec2[4] =  7;

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      //checkCapacity( vec2, 5UL );
      //checkNonZeros( vec2, 2UL );

      if( vec2[0] != 0 || vec2[1] != 0 || vec2[2] != -4 || vec2[3] != -18 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 0 -4 -18 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all CustomVector (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the CustomVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testScaling()
{
   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "CustomVector self-scaling (v*=s)";

      VT vec( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec[0] =  1;
      vec[1] =  0;
      vec[2] = -2;
      vec[3] =  3;
      vec[4] =  0;

      vec *= 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 0 || vec[2] != -4 || vec[3] != 6 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 0 -4 6 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "CustomVector self-scaling (v=v*s)";

      VT vec( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec[0] =  1;
      vec[1] =  0;
      vec[2] = -2;
      vec[3] =  3;
      vec[4] =  0;

      vec = vec * 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 0 || vec[2] != -4 || vec[3] != 6 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 0 -4 6 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "CustomVector self-scaling (v=s*v)";

      VT vec( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec[0] =  1;
      vec[1] =  0;
      vec[2] = -2;
      vec[3] =  3;
      vec[4] =  0;

      vec = 2 * vec;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 0 || vec[2] != -4 || vec[3] != 6 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 0 -4 6 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "CustomVector self-scaling (v/=s)";

      VT vec( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec[0] =  2;
      vec[1] =  0;
      vec[2] = -4;
      vec[3] =  6;
      vec[4] =  0;

      vec /= 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 0 || vec[2] != -2 || vec[3] != 3 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 -2 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "CustomVector self-scaling (v=v/s)";

      VT vec( blaze::allocate<int>( 5UL ), 5UL, blaze::Deallocate() );
      vec[0] =  2;
      vec[1] =  0;
      vec[2] = -4;
      vec[3] =  6;
      vec[4] =  0;

      vec = vec / 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 0 || vec[2] != -2 || vec[3] != 3 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 -2 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // CustomVector::scale()
   //=====================================================================================

   {
      test_ = "CustomVector::scale() (int)";

      // Initialization check
      VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 3;
      vec[3] = 4;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the vector
      vec.scale( 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 4 || vec[2] != 6 || vec[3] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 4 6 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the vector
      vec.scale( 0.5 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CustomVector::scale() (complex)";

      using blaze::complex;
      using blaze::aligned;
      using blaze::unpadded;
      using blaze::rowVector;

      typedef blaze::CustomVector<complex<float>,aligned,unpadded,rowVector>  AlignedUnpadded;
      AlignedUnpadded vec( blaze::allocate< complex<float> >( 2UL ), 2UL, blaze::Deallocate() );
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
/*!\brief Test of the CustomVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the CustomVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void AlignedUnpaddedTest::testSubscript()
{
   test_ = "CustomVector::operator[]";

   // Assignment to the element at index 2
   VT vec( blaze::allocate<int>( 7UL ), 7UL, blaze::Deallocate() );
   reset( vec );
   vec[2] = 1;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
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

   // Assignment to the element at index 5
   vec[5] = 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
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

   // Assignment to the element at index 3
   vec[3] = 3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
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

   // Assignment to the element at index 0
   vec[0] = 4;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
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

   // Addition assignment to the element at index 2
   vec[2] += vec[3];

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 4UL );

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
   checkNonZeros( vec, 5UL );

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
   checkNonZeros( vec, 5UL );

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
   checkNonZeros( vec, 5UL );

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
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the CustomVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the CustomVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void AlignedUnpaddedTest::testAt()
{
   test_ = "CustomVector::at()";

   // Assignment to the element at index 2
   VT vec( blaze::allocate<int>( 7UL ), 7UL, blaze::Deallocate() );
   reset( vec );
   vec.at(2) = 1;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 1UL );

   if( vec.at(2) != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 5
   vec.at(5) = 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 2UL );

   if( vec.at(2) != 1 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 3
   vec.at(3) = 3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 3UL );

   if( vec.at(2) != 1 || vec.at(3) != 3 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 0
   vec.at(0) = 4;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 4UL );

   if( vec.at(0) != 4 || vec.at(2) != 1 || vec.at(3) != 3 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Addition assignment to the element at index 2
   vec.at(2) += vec.at(3);

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 4UL );

   if( vec.at(0) != 4 || vec.at(2) != 4 || vec.at(3) != 3 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 4 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Subtraction assignment to the element at index 1
   vec.at(1) -= vec.at(5);

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 4 || vec.at(3) != 3 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 3
   vec.at(3) *= -3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 4 || vec.at(3) != -9 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 -9 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Division assignment to the element at index 2
   vec.at(2) /= 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 2 || vec.at(3) != -9 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Attempt to assign to the element at index 7
   try {
      vec.at(7) = 2;

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomVector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the CustomVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testIterator()
{
   typedef VT::Iterator       Iterator;
   typedef VT::ConstIterator  ConstIterator;

   VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
   vec[0] =  1;
   vec[1] =  0;
   vec[2] = -2;
   vec[3] = -3;

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      Iterator it = Iterator();

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

      ConstIterator it = ConstIterator();

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
/*!\brief Test of the \c nonZeros() member function of the CustomVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the CustomVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testNonZeros()
{
   test_ = "CustomVector::nonZeros()";

   {
      VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
      reset( vec );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 0;
      vec[3] = 3;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the CustomVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the CustomVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testReset()
{
   using blaze::reset;


   //=====================================================================================
   // CustomVector::reset()
   //=====================================================================================

   {
      test_ = "CustomVector::reset()";

      // Initialization check
      VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 3;
      vec[3] = 4;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( vec[2] );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the vector
      reset( vec );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // CustomVector::reset( Type*, size_t, size_t )
   //=====================================================================================

   {
      test_ = "CustomVector::reset( Type*, size_t, size_t )";

      VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 3;
      vec[3] = 4;

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 27UL ) );
      vec.reset( array.get(), 27UL );

      checkSize    ( vec, 27UL );
      checkCapacity( vec, 27UL );
   }


   //=====================================================================================
   // CustomVector::reset( Type*, size_t, size_t, Deleter )
   //=====================================================================================

   {
      test_ = "CustomVector::reset( Type*, size_t, size_t, Deleter )";

      VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 3;
      vec[3] = 4;

      vec.reset( blaze::allocate<int>( 27UL ), 27UL, blaze::Deallocate() );

      checkSize    ( vec, 27UL );
      checkCapacity( vec, 27UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the CustomVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the CustomVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testClear()
{
   using blaze::clear;

   test_ = "CustomVector::clear()";

   // Initialization check
   VT vec( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
   vec[0] = 1;
   vec[1] = 2;
   vec[2] = 3;
   vec[3] = 4;

   checkSize    ( vec, 4UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 3 4 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Clearing a single element
   clear( vec[2] );

   checkSize    ( vec, 4UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 3UL );

   if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Clear operation failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 0 4 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Clearing the vector
   clear( vec );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the CustomVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the CustomVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testSwap()
{
   test_ = "CustomVector swap";

   VT vec1( blaze::allocate<int>( 3UL ), 3UL, blaze::Deallocate() );
   vec1[0] = 1;
   vec1[1] = 2;
   vec1[2] = 3;

   VT vec2( blaze::allocate<int>( 4UL ), 4UL, blaze::Deallocate() );
   vec2[0] = 4;
   vec2[1] = 3;
   vec2[2] = 2;
   vec2[3] = 1;

   swap( vec1, vec2 );

   checkSize    ( vec1, 4UL );
   checkCapacity( vec1, 4UL );
   checkNonZeros( vec1, 4UL );

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
   checkNonZeros( vec2, 3UL );

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
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the CustomVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the CustomVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testIsDefault()
{
   using blaze::isDefault;

   test_ = "isDefault() function";

   // isDefault with vector of size 0
   {
      VT vec;

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
      VT vec( blaze::allocate<int>( 3UL ), 3UL, blaze::Deallocate() );
      reset( vec );

      if( isDefault( vec[1] ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector element: " << vec[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( vec ) != false ) {
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
      VT vec( blaze::allocate<int>( 3UL ), 3UL, blaze::Deallocate() );
      reset( vec );
      vec[1] = 1;

      if( isDefault( vec[1] ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector element: " << vec[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

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

} // namespace customvector

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
   std::cout << "   Running aligned/unpadded CustomVector class test..." << std::endl;

   try
   {
      RUN_CUSTOMVECTOR_ALIGNED_UNPADDED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aligned/unpadded CustomVector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
