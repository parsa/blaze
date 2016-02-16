//=================================================================================================
/*!
//  \file src/mathtest/custommatrix/AlignedPaddedTest.cpp
//  \brief Source file for the aligned/padded CustomMatrix class test
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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/Functions.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/policies/ArrayDelete.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blaze/util/typetraits/IsVectorizable.h>
#include <blaze/util/UniqueArray.h>
#include <blazetest/mathtest/custommatrix/AlignedPaddedTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace mathtest {

namespace custommatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the CustomMatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
AlignedPaddedTest::AlignedPaddedTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testScaling();
   testFunctionCall();
   testAt();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testTranspose();
   testCTranspose();
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
/*!\brief Test of the CustomMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the CustomMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix default constructor";

      MT mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Row-major constructor ( Type*, size_t, size_t, size_t )
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix constructor ( Type*, size_t, size_t, size_t )";

      // Constructor a 2x3 custom matrix
      {
         blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 32UL ) );
         MT mat( array.get(), 2UL, 3UL, 16UL );

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 32UL );
      }

      // Trying to construct a custom matrix with invalid array of elements
      try {
         MT mat( NULL, 0UL, 0UL, 0UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom matrix with a NULL pointer succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom matrix with invalid alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 33UL ) );
            MT mat( array.get()+1UL, 2UL, 2UL, 16UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid row alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 30UL ) );
            MT mat( array.get(), 2UL, 2UL, 15UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid row alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid padding
      if( blaze::IsVectorizable<int>::value )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 6UL ) );
            MT mat( array.get(), 2UL, 2UL, 3UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid padding succeeded\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Row-major constructor ( Type*, size_t, size_t, size_t, Deleter )
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix constructor ( Type*, size_t, size_t, size_t, Deleter )";

      // Constructor a 2x3 custom matrix
      {
         MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 32UL );
      }

      // Trying to construct a custom matrix with invalid array of elements
      try {
         MT mat( NULL, 0UL, 0UL, 0UL, blaze::Deallocate() );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom matrix with a NULL pointer succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom matrix with invalid alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 33UL ) );
            MT mat( array.get()+1UL, 2UL, 2UL, 16UL, blaze::Deallocate() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid row alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 30UL ) );
            MT mat( array.get(), 2UL, 2UL, 15UL, blaze::Deallocate() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid row alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid padding
      if( blaze::IsVectorizable<int>::value )
      {
         try {
            MT mat( blaze::allocate<int>( 6UL ), 2UL, 2UL, 3UL, blaze::Deallocate() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid padding succeeded\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix copy constructor (0x0)";

      MT mat1;
      MT mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix copy constructor (0x3)";

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 16UL ) );
      MT mat1( array.get(), 0UL, 3UL, 16UL );
      MT mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix copy constructor (2x0)";

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 10UL ) );
      MT mat1( array.get(), 2UL, 0UL, 0UL );
      MT mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix copy constructor (2x3)";

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 32UL ) );
      MT mat1( array.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      MT mat2( mat1 );

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix default constructor";

      OMT mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major constructor ( Type*, size_t, size_t, size_t )
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix constructor ( Type*, size_t, size_t, size_t )";

      // Constructor a 2x3 custom matrix
      {
         blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 48UL ) );
         OMT mat( array.get(), 2UL, 3UL, 16UL );

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
      }

      // Trying to construct a custom matrix with invalid array of elements
      try {
         OMT mat( NULL, 0UL, 0UL, 0UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom matrix with a NULL pointer succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom matrix with invalid alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 33UL ) );
            OMT mat( array.get()+1UL, 2UL, 2UL, 16UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid column alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 30UL ) );
            OMT mat( array.get(), 2UL, 2UL, 15UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid column alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid padding
      if( blaze::IsVectorizable<int>::value )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 6UL ) );
            OMT mat( array.get(), 2UL, 2UL, 3UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid padding succeeded\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major constructor ( Type*, size_t, size_t, size_t, Deleter )
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix constructor ( Type*, size_t, size_t, size_t, Deleter )";

      // Constructor a 2x3 custom matrix
      {
         OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
      }

      // Trying to construct a custom matrix with invalid array of elements
      try {
         OMT mat( NULL, 0UL, 0UL, 0UL, blaze::Deallocate() );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom matrix with a NULL pointer succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom matrix with invalid alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 33UL ) );
            OMT mat( array.get()+1UL, 2UL, 2UL, 16UL, blaze::Deallocate() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid column alignment
      if( blaze::AlignmentOf<int>::value > 1UL )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 30UL ) );
            OMT mat( array.get(), 2UL, 2UL, 15UL, blaze::Deallocate() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid column alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid padding
      if( blaze::IsVectorizable<int>::value )
      {
         try {
            blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 6UL ) );
            OMT mat( array.get(), 2UL, 2UL, 3UL, blaze::Deallocate() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid padding succeeded\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix copy constructor (0x0)";

      OMT mat1;
      OMT mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix copy constructor (0x3)";

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 10UL ) );
      OMT mat1( array.get(), 0UL, 3UL, 0UL );
      OMT mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix copy constructor (2x0)";

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 16UL ) );
      OMT mat1( array.get(), 2UL, 0UL, 16UL );
      OMT mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix copy constructor (2x3)";

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 48UL ) );
      OMT mat1( array.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      OMT mat2( mat1 );

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the CustomMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix homogeneous assignment";

      MT mat( blaze::allocate<int>( 48UL ), 3UL, 4UL, 16UL, blaze::Deallocate() );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 4UL );
      checkNonZeros( mat,  1UL, 4UL );
      checkNonZeros( mat,  2UL, 4UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array assignment
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix array assignment";

      const int array[2UL][3UL] = { { 1, 2, 3 }, { 4, 5, 6 } };
      MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = array;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix copy assignment";

      MT mat1( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<unsigned int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<unsigned int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1(0,0) = 1U;
      mat1(0,1) = 2U;
      mat1(0,2) = 3U;
      mat1(1,0) = 4U;
      mat1(1,1) = 5U;
      mat1(1,2) = 6U;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment stress test (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      const unsigned int min( randmin );
      const unsigned int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         typedef blaze::CustomMatrix<unsigned int,aligned,padded,rowMajor>  AlignedPadded;
         AlignedPadded mat1( blaze::allocate<unsigned int>( rows*spacing ), rows, columns, spacing, blaze::Deallocate() );
         randomize( mat1, min, max );

         MT mat2( blaze::allocate<int>( rows*spacing ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment stress test (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
         blaze::UniqueArray<int> array( new int[rows*columns+1UL] );
         UnalignedUnpadded mat1( array.get()+1UL, rows, columns );
         randomize( mat1, min, max );

         MT mat2( blaze::allocate<int>( rows*spacing ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment stress test (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows    ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing1( blaze::nextMultiple<size_t>( rows, 16UL ) );
         const size_t spacing2( blaze::nextMultiple<size_t>( columns, 16UL ) );

         typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
         AlignedPadded mat1( blaze::allocate<int>( spacing1*columns ), rows, columns, spacing1, blaze::Deallocate() );
         randomize( mat1, min, max );

         MT mat2( blaze::allocate<int>( rows*spacing2 ), rows, columns, spacing2, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment stress test (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
         blaze::UniqueArray<int> array( new int[rows*columns+1UL] );
         UnalignedUnpadded mat1( array.get()+1UL, rows, columns );
         randomize( mat1, min, max );

         MT mat2( blaze::allocate<int>( rows*spacing ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment stress test";

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         blaze::CompressedMatrix<int,blaze::rowMajor> mat1( rows, columns );
         randomize( mat1, min, max );

         MT mat2( blaze::allocate<int>( rows*spacing ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment stress test";

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         blaze::CompressedMatrix<int,blaze::columnMajor> mat1( rows, columns );
         randomize( mat1, min, max );

         MT mat2( blaze::allocate<int>( rows*spacing ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix homogeneous assigment";

      OMT mat( blaze::allocate<int>( 64UL ), 3UL, 4UL, 16UL, blaze::Deallocate() );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 64UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );
      checkNonZeros( mat,  2UL, 3UL );
      checkNonZeros( mat,  3UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array assignment
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix array initialization constructor";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = array;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix copy assignment";

      OMT mat1( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment stress test (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows    ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing1( blaze::nextMultiple<size_t>( columns, 16UL ) );
         const size_t spacing2( blaze::nextMultiple<size_t>( rows, 16UL ) );

         typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
         AlignedPadded mat1( blaze::allocate<int>( rows*spacing1 ), rows, columns, spacing1, blaze::Deallocate() );
         randomize( mat1, min, max );

         OMT mat2( blaze::allocate<int>( spacing2*columns ), rows, columns, spacing2, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment stress test (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
         blaze::UniqueArray<int> array( new int[rows*columns+1UL] );
         UnalignedUnpadded mat1( array.get()+1UL, rows, columns );
         randomize( mat1, min, max );

         OMT mat2( blaze::allocate<int>( spacing*columns ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<unsigned int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<unsigned int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1(0,0) = 1U;
      mat1(0,1) = 2U;
      mat1(0,2) = 3U;
      mat1(1,0) = 4U;
      mat1(1,1) = 5U;
      mat1(1,2) = 6U;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment stress test (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      const unsigned int min( randmin );
      const unsigned int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         typedef blaze::CustomMatrix<unsigned int,aligned,padded,columnMajor>  AlignedPadded;
         AlignedPadded mat1( blaze::allocate<unsigned int>( spacing*columns ), rows, columns, spacing, blaze::Deallocate() );
         randomize( mat1, min, max );

         OMT mat2( blaze::allocate<int>( spacing*columns ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment stress test (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
         blaze::UniqueArray<int> array( new int[rows*columns+1UL] );
         UnalignedUnpadded mat1( array.get()+1UL, rows, columns );
         randomize( mat1, min, max );

         OMT mat2( blaze::allocate<int>( spacing*columns ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 1UL );
      checkNonZeros( mat2,  2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment stress test";

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         blaze::CompressedMatrix<int,blaze::rowMajor> mat1( rows, columns );
         randomize( mat1, min, max );

         OMT mat2( blaze::allocate<int>( spacing*columns ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 1UL );
      checkNonZeros( mat2,  2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment stress test";

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         blaze::CompressedMatrix<int,blaze::columnMajor> mat1( rows, columns );
         randomize( mat1, min, max );

         OMT mat2( blaze::allocate<int>( spacing*columns ), rows, columns, spacing, blaze::Deallocate() );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the CustomMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the CustomMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      MT mat2( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[7UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      OMT mat2( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the CustomMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[10UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 3UL, 3UL );
      mat1 = 0;
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[10UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 3UL, 3UL );
      mat1 = 0;
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 3UL, 3UL, 5UL );
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 3UL, 3UL, 5UL );
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      MT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[10UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 3UL, 3UL );
      mat1 = 0;
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat1( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat1 = 0;
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      blaze::UniqueArray<int> array( new int[10UL] );
      UnalignedUnpadded mat1( array.get()+1UL, 3UL, 3UL );
      mat1 = 0;
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 3UL, 3UL, 5UL );
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 3UL, 3UL, 5UL );
      mat1(0,1) = 2;
      mat1(1,0) = 1;
      mat1(1,1) = 3;
      mat1(1,2) = 4;
      mat1(2,2) = 5;

      OMT mat2( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat2 = 0;
      mat2(0,0) = 1;
      mat2(0,2) = 2;
      mat2(1,1) = 3;
      mat2(2,0) = 4;
      mat2(2,2) = 5;

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 10 )\n( 3 9 12 )\n( 0 8 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all CustomMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the CustomMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      MT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  1;
      mat(2,0) = -2;
      mat(2,2) =  3;

      mat *= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      MT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  1;
      mat(2,0) = -2;
      mat(2,2) =  3;

      mat = mat * 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      MT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  1;
      mat(2,0) = -2;
      mat(2,2) =  3;

      mat = 2 * mat;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      MT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  2;
      mat(2,0) = -4;
      mat(2,2) =  6;

      mat /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 1 ||
          mat(2,0) != -2 || mat(2,1) != 0 || mat(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      MT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  2;
      mat(2,0) = -4;
      mat(2,2) =  6;

      mat = mat / 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 1 ||
          mat(2,0) != -2 || mat(2,1) != 0 || mat(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major CustomMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix::scale() (int)";

      // Initialization check
      MT mat( blaze::allocate<int>( 48UL ), 3UL, 2UL, 16UL, blaze::Deallocate() );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(1,0) = 3;
      mat(1,1) = 4;
      mat(2,0) = 5;
      mat(2,1) = 6;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 ||
          mat(1,0) != 3 || mat(1,1) != 4 ||
          mat(2,0) != 5 || mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 )\n( 3 4 )\n( 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      mat.scale( 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,0) !=  2 || mat(0,1) !=  4 ||
          mat(1,0) !=  6 || mat(1,1) !=  8 ||
          mat(2,0) != 10 || mat(2,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  2  4 )\n(  6  8 )\n( 10 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      mat.scale( 0.5 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 ||
          mat(1,0) != 3 || mat(1,1) != 4 ||
          mat(2,0) != 5 || mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 )\n( 3 4 )\n( 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major CustomMatrix::scale() (complex)";

      using blaze::complex;
      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<complex<float>,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate< complex<float> >( 32UL ), 2UL, 2UL, 16UL, blaze::Deallocate() );
      mat(0,0) = complex<float>( 1.0F, 0.0F );
      mat(0,1) = complex<float>( 2.0F, 0.0F );
      mat(1,0) = complex<float>( 3.0F, 0.0F );
      mat(1,1) = complex<float>( 4.0F, 0.0F );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );

      if( mat(0,0) != complex<float>( 3.0F, 0.0F ) || mat(0,1) != complex<float>(  6.0F, 0.0F ) ||
          mat(1,0) != complex<float>( 9.0F, 0.0F ) || mat(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( ( 3,0) ( 6,0)\n( 9,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      OMT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  1;
      mat(2,0) = -2;
      mat(2,2) =  3;

      mat *= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      OMT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  1;
      mat(2,0) = -2;
      mat(2,2) =  3;

      mat = mat * 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      OMT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  1;
      mat(2,0) = -2;
      mat(2,2) =  3;

      mat = 2 * mat;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      OMT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  2;
      mat(2,0) = -4;
      mat(2,2) =  6;

      mat /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 1 ||
          mat(2,0) != -2 || mat(2,1) != 0 || mat(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      OMT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,2) =  2;
      mat(2,0) = -4;
      mat(2,2) =  6;

      mat = mat / 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 1 ||
          mat(2,0) != -2 || mat(2,1) != 0 || mat(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major CustomMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix::scale() (int)";

      // Initialization check
      OMT mat( blaze::allocate<int>( 32UL ), 3UL, 2UL, 16UL, blaze::Deallocate() );
      mat(0,0) = 1;
      mat(0,1) = 4;
      mat(1,0) = 2;
      mat(1,1) = 5;
      mat(2,0) = 3;
      mat(2,1) = 6;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 4 ||
          mat(1,0) != 2 || mat(1,1) != 5 ||
          mat(2,0) != 3 || mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 4 )\n( 2 5 )\n( 3 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      mat.scale( 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) !=  8 ||
          mat(1,0) != 4 || mat(1,1) != 10 ||
          mat(2,0) != 6 || mat(2,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  2  8 )\n(  4 10 )\n(  6 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      mat.scale( 0.5 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 4 ||
          mat(1,0) != 2 || mat(1,1) != 5 ||
          mat(2,0) != 3 || mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 4 )\n( 2 5 )\n( 3 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major CustomMatrix::scale() (complex)";

      using blaze::complex;
      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<complex<float>,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate< complex<float> >( 32UL ), 2UL, 2UL, 16UL, blaze::Deallocate() );
      mat(0,0) = complex<float>( 1.0F, 0.0F );
      mat(0,1) = complex<float>( 2.0F, 0.0F );
      mat(1,0) = complex<float>( 3.0F, 0.0F );
      mat(1,1) = complex<float>( 4.0F, 0.0F );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );

      if( mat(0,0) != complex<float>( 3.0F, 0.0F ) || mat(0,1) != complex<float>(  6.0F, 0.0F ) ||
          mat(1,0) != complex<float>( 9.0F, 0.0F ) || mat(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( ( 3,0) ( 6,0)\n( 9,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the CustomMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void AlignedPaddedTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix::operator()";

      // Assignment to the element (2,1)
      MT mat( blaze::allocate<int>( 48UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(2,1) = 1;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  1UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 0UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (1,4)
      mat(1,4) = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  2UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat(1,4) != 2 || mat(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (0,3)
      mat(0,3) = 3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  3UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat(0,3) != 3 || mat(1,4) != 2 || mat(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (2,2)
      mat(2,2) = 4;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,3) != 3 || mat(1,4) != 2 || mat(2,1) != 1 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element (2,1)
      mat(2,1) += mat(0,3);

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,3) != 3 || mat(1,4) != 2 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element (1,0)
      mat(1,0) -= mat(1,4);

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,3) != 3 || mat(1,0) != -2 || mat(1,4) != 2 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 3 0 )\n( -2 0 0 0 2 )\n(  0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element (0,3)
      mat(0,3) *= -3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,3) != -9 || mat(1,0) != -2 || mat(1,4) != 2 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 4 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element (2,1)
      mat(2,1) /= 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,3) != -9 || mat(1,0) != -2 || mat(1,4) != 2 || mat(2,1) != 2 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix::operator()";

      // Assignment to the element (2,1)
      OMT mat( blaze::allocate<int>( 80UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(2,1) = 1;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  1UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 0UL );
      checkNonZeros( mat,  3UL, 0UL );
      checkNonZeros( mat,  4UL, 0UL );

      if( mat(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (1,4)
      mat(1,4) = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  2UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 0UL );
      checkNonZeros( mat,  3UL, 0UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat(2,1) != 1 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (0,3)
      mat(0,3) = 3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  3UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 0UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat(2,1) != 1 || mat(1,4) != 2 || mat(0,3) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (2,2)
      mat(2,2) = 4;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat(2,1) != 1 || mat(1,4) != 2 || mat(0,3) != 3 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element (2,1)
      mat(2,1) += mat(0,3);

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat(2,1) != 4 || mat(2,2) != 4 || mat(0,3) != 3 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element (1,0)
      mat(1,0) -= mat(1,4);

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat(1,0) != -2 || mat(2,1) != 4 || mat(2,2) != 4 || mat(0,3) != 3 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 3 0 )\n( -2 0 0 0 2 )\n(  0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element (0,3)
      mat(0,3) *= -3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat(1,0) != -2 || mat(2,1) != 4 || mat(2,2) != 4 || mat(0,3) != -9 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 4 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element (2,1)
      mat(2,1) /= 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat(1,0) != -2 || mat(2,1) != 2 || mat(2,2) != 4 || mat(0,3) != -9 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the CustomMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the CustomMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void AlignedPaddedTest::testAt()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix::at()";

      // Assignment to the element (2,1)
      MT mat( blaze::allocate<int>( 48UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat.at(2,1) = 1;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  1UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 0UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat.at(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (1,4)
      mat.at(1,4) = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  2UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat.at(1,4) != 2 || mat.at(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (0,3)
      mat.at(0,3) = 3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  3UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat.at(0,3) != 3 || mat.at(1,4) != 2 || mat.at(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (2,2)
      mat.at(2,2) = 4;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat.at(0,3) != 3 || mat.at(1,4) != 2 || mat.at(2,1) != 1 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element (2,1)
      mat.at(2,1) += mat.at(0,3);

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat.at(0,3) != 3 || mat.at(1,4) != 2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element (1,0)
      mat.at(1,0) -= mat.at(1,4);

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat.at(0,3) != 3 || mat.at(1,0) != -2 || mat.at(1,4) != 2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 3 0 )\n( -2 0 0 0 2 )\n(  0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element (0,3)
      mat.at(0,3) *= -3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat.at(0,3) != -9 || mat.at(1,0) != -2 || mat.at(1,4) != 2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 4 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element (2,1)
      mat.at(2,1) /= 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat.at(0,3) != -9 || mat.at(1,0) != -2 || mat.at(1,4) != 2 || mat.at(2,1) != 2 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Attempt to assign to the element (3,0)
      try {
         mat.at(3,0) = 2;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}

      // Attempt to assign to the element (0,5)
      try {
         mat.at(0,5) = 2;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix::at()";

      // Assignment to the element (2,1)
      OMT mat( blaze::allocate<int>( 80UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat.at(2,1) = 1;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  1UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 0UL );
      checkNonZeros( mat,  3UL, 0UL );
      checkNonZeros( mat,  4UL, 0UL );

      if( mat.at(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (1,4)
      mat.at(1,4) = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  2UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 0UL );
      checkNonZeros( mat,  3UL, 0UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat.at(2,1) != 1 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (0,3)
      mat.at(0,3) = 3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  3UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 0UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat.at(2,1) != 1 || mat.at(1,4) != 2 || mat.at(0,3) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (2,2)
      mat.at(2,2) = 4;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat.at(2,1) != 1 || mat.at(1,4) != 2 || mat.at(0,3) != 3 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element (2,1)
      mat.at(2,1) += mat.at(0,3);

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat.at(2,1) != 4 || mat.at(2,2) != 4 || mat.at(0,3) != 3 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element (1,0)
      mat.at(1,0) -= mat.at(1,4);

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat.at(1,0) != -2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 || mat.at(0,3) != 3 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 3 0 )\n( -2 0 0 0 2 )\n(  0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element (0,3)
      mat.at(0,3) *= -3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat.at(1,0) != -2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 || mat.at(0,3) != -9 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 4 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element (2,1)
      mat.at(2,1) /= 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 1UL );

      if( mat.at(1,0) != -2 || mat.at(2,1) != 2 || mat.at(2,2) != 4 || mat.at(0,3) != -9 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Attempt to assign to the element (3,0)
      try {
         mat.at(3,0) = 2;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}

      // Attempt to assign to the element (0,5)
      try {
         mat.at(0,5) = 2;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the CustomMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      typedef MT::Iterator       Iterator;
      typedef MT::ConstIterator  ConstIterator;

      MT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(0,1) =  1;
      mat(1,0) = -2;
      mat(1,2) = -3;
      mat(2,1) =  4;
      mat(2,2) =  5;

      // Testing the Iterator default constructor
      {
         test_ = "Row-major Iterator default constructor";

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
         test_ = "Row-major ConstIterator default constructor";

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
         test_ = "Row-major Iterator/ConstIterator conversion";

         ConstIterator it( begin( mat, 1UL ) );

         if( it == end( mat, 1UL ) || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator
      {
         test_ = "Row-major Iterator subtraction";

         const size_t number( end( mat, 0UL ) - begin( mat, 0UL ) );

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

      // Counting the number of elements in 1st row via ConstIterator
      {
         test_ = "Row-major ConstIterator subtraction";

         const size_t number( cend( mat, 1UL ) - cbegin( mat, 1UL ) );

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
         test_ = "Row-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( mat, 2UL ) );
         ConstIterator end( cend( mat, 2UL ) );

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 2UL;

         if( it == end || *it != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 2UL;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar subtraction failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = 3UL + it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Scalar/iterator addition failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Row-major assignment via Iterator";

         int value = 7;

         for( Iterator it=begin( mat, 2UL ); it!=end( mat, 2UL ); ++it ) {
            *it = value++;
         }

         if( mat(0,0) !=  0 || mat(0,1) != 1 || mat(0,2) !=  0 ||
             mat(1,0) != -2 || mat(1,1) != 0 || mat(1,2) != -3 ||
             mat(2,0) !=  7 || mat(2,1) != 8 || mat(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it += value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 0 ||
             mat(1,0) != 2 || mat(1,1) != 5 || mat(1,2) != 3 ||
             mat(2,0) != 7 || mat(2,1) != 8 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 2 5 3 )\n( 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it -= value++;
         }

         if( mat(0,0) !=  0 || mat(0,1) != 1 || mat(0,2) !=  0 ||
             mat(1,0) != -2 || mat(1,1) != 0 || mat(1,2) != -3 ||
             mat(2,0) !=  7 || mat(2,1) != 8 || mat(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         int value = 2;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it *= value++;
         }

         if( mat(0,0) !=  0 || mat(0,1) != 1 || mat(0,2) !=   0 ||
             mat(1,0) != -4 || mat(1,1) != 0 || mat(1,2) != -12 ||
             mat(2,0) !=  7 || mat(2,1) != 8 || mat(2,2) !=   9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n(  0  1   0 )\n( -4  0 -12 )\n(  7  8   9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it /= 2;
         }

         if( mat(0,0) !=  0 || mat(0,1) != 1 || mat(0,2) !=  0 ||
             mat(1,0) != -2 || mat(1,1) != 0 || mat(1,2) != -6 ||
             mat(2,0) !=  7 || mat(2,1) != 8 || mat(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -6 )\n(  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      typedef OMT::Iterator       Iterator;
      typedef OMT::ConstIterator  ConstIterator;


      OMT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,0) =  1;
      mat(0,1) = -2;
      mat(2,1) = -3;
      mat(1,2) =  4;
      mat(2,2) =  5;

      // Testing the Iterator default constructor
      {
         test_ = "Column-major Iterator default constructor";

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
         test_ = "Column-major ConstIterator default constructor";

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
         test_ = "Column-major Iterator/ConstIterator conversion";

         ConstIterator it( begin( mat, 1UL ) );

         if( it == end( mat, 1UL ) || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator
      {
         test_ = "Column-major Iterator subtraction";

         const size_t number( end( mat, 0UL ) - begin( mat, 0UL ) );

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

      // Counting the number of elements in 1st row via ConstIterator
      {
         test_ = "Column-major ConstIterator subtraction";

         const size_t number( cend( mat, 1UL ) - cbegin( mat, 1UL ) );

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
         test_ = "Column-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( mat, 2UL ) );
         ConstIterator end( cend( mat, 2UL ) );

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 2UL;

         if( it == end || *it != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 2UL;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar subtraction failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = 3UL + it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Scalar/iterator addition failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Column-major assignment via Iterator";

         int value = 7;

         for( Iterator it=begin( mat, 2UL ); it!=end( mat, 2UL ); ++it ) {
            *it = value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != -2 || mat(0,2) != 7 ||
             mat(1,0) != 1 || mat(1,1) !=  0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != -3 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 -2  7 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it += value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != 2 || mat(0,2) != 7 ||
             mat(1,0) != 1 || mat(1,1) != 5 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != 3 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 2 7 )\n( 1 5 8 )\n( 0 3 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it -= value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != -2 || mat(0,2) != 7 ||
             mat(1,0) != 1 || mat(1,1) !=  0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != -3 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 -2  7 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         int value = 2;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it *= value++;
         }

         if( mat(0,0) != 0 || mat(0,1) !=  -4 || mat(0,2) != 7 ||
             mat(1,0) != 1 || mat(1,1) !=   0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != -12 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 -2  7 )\n( 1  0  8 )\n( 0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it /= 2;
         }

         if( mat(0,0) != 0 || mat(0,1) != -2 || mat(0,2) != 7 ||
             mat(1,0) != 1 || mat(1,1) !=  0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != -6 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 -2  7 )\n( 1  0  8 )\n( 0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the CustomMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the CustomMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix::nonZeros()";

      {
         MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
         mat = 0;

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 32UL );
         checkNonZeros( mat,  0UL );
         checkNonZeros( mat,  0UL, 0UL );
         checkNonZeros( mat,  1UL, 0UL );

         if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
             mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
         mat = 0;
         mat(0,1) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 32UL );
         checkNonZeros( mat,  3UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );

         if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 2 ||
             mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 1 2 )\n( 0 3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix::nonZeros()";

      {
         OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
         mat = 0;

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  0UL );
         checkNonZeros( mat,  0UL, 0UL );
         checkNonZeros( mat,  1UL, 0UL );
         checkNonZeros( mat,  2UL, 0UL );

         if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
             mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
         mat = 0;
         mat(0,1) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  3UL );
         checkNonZeros( mat,  0UL, 0UL );
         checkNonZeros( mat,  1UL, 2UL );
         checkNonZeros( mat,  2UL, 1UL );

         if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 2 ||
             mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 1 2 )\n( 0 3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the CustomMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the CustomMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testReset()
{
   using blaze::reset;


   //=====================================================================================
   // Row-major CustomMatrix::reset()
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix::reset()";

      // Initialization check
      MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(0,2) = 3;
      mat(1,0) = 4;
      mat(1,1) = 5;
      mat(1,2) = 6;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( mat(0,2) );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( mat, 1UL );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  2UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 0UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( mat );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  0UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 0UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major CustomMatrix::reset( Type*, size_t, size_t, size_t )
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix::reset( Type*, size_t, size_t, size_t )";

      MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 2;

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 48UL ) );
      mat.reset( array.get(), 3UL, 5UL, 16UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
   }


   //=====================================================================================
   // Row-major CustomMatrix::reset( Type*, size_t, size_t, size_t, Deleter )
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix::reset( Type*, size_t, size_t, size_t, Deleter )";

      MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 2;

      mat.reset( blaze::allocate<int>( 48UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 48UL );
   }


   //=====================================================================================
   // Column-major CustomMatrix::reset()
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix::reset()";

      // Initialization check
      OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(0,2) = 3;
      mat(1,0) = 4;
      mat(1,1) = 5;
      mat(1,2) = 6;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( mat(0,2) );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting column 1
      reset( mat, 1UL );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  3UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 0UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( mat );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  0UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 0UL );
      checkNonZeros( mat,  2UL, 0UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major CustomMatrix::reset( Type*, size_t, size_t, size_t )
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix::reset( Type*, size_t, size_t, size_t )";

      OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 2;

      blaze::UniqueArray<int,blaze::Deallocate> array( blaze::allocate<int>( 80UL ) );
      mat.reset( array.get(), 3UL, 5UL, 16UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
   }


   //=====================================================================================
   // Column-major CustomMatrix::reset( Type*, size_t, size_t, size_t, Deleter )
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix::reset( Type*, size_t, size_t, size_t, Deleter )";

      OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 2;

      mat.reset( blaze::allocate<int>( 80UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 80UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the CustomMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the CustomMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testClear()
{
   using blaze::clear;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix::clear()";

      // Initialization check
      MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(0,2) = 3;
      mat(1,0) = 4;
      mat(1,1) = 5;
      mat(1,2) = 6;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( mat(0,2) );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( mat );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix::clear()";

      // Initialization check
      OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(0,2) = 3;
      mat(1,0) = 4;
      mat(1,1) = 5;
      mat(1,2) = 6;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( mat(0,2) );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  5UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( mat );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the CustomMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the CustomMatrix
// class template. Additionally, it performs a test of self-transpose via the \c trans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via transpose()";

      // Self-transpose of a 3x3 matrix
      {
         MT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
         mat(0,0) = 1;
         mat(0,1) = 0;
         mat(0,2) = 2;
         mat(1,0) = 0;
         mat(1,1) = 3;
         mat(1,2) = 0;
         mat(2,0) = 4;
         mat(2,1) = 0;
         mat(2,2) = 5;

         transpose( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  5UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );
         checkNonZeros( mat,  2UL, 2UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 4 ||
             mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 ||
             mat(2,0) != 2 || mat(2,1) != 0 || mat(2,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 4 )\n( 0 3 0 )\n( 2 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Try to self-transpose a 3x5 matrix
      try {
         MT mat( blaze::allocate<int>( 48UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );

         transpose( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Self-transpose of a non-square matrix succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::logic_error& ) {}
   }

   {
      test_ = "Row-major self-transpose via trans()";

      // Self-transpose of a 3x3 matrix
      {
         MT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
         mat(0,0) = 1;
         mat(0,1) = 0;
         mat(0,2) = 2;
         mat(1,0) = 0;
         mat(1,1) = 3;
         mat(1,2) = 0;
         mat(2,0) = 4;
         mat(2,1) = 0;
         mat(2,2) = 5;

         mat = trans( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  5UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );
         checkNonZeros( mat,  2UL, 2UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 4 ||
             mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 ||
             mat(2,0) != 2 || mat(2,1) != 0 || mat(2,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 4 )\n( 0 3 0 )\n( 2 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Try to self-transpose a 3x5 matrix
      try {
         MT mat( blaze::allocate<int>( 48UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );

         mat = trans( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Self-transpose of a non-square matrix succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      // Self-transpose of a 3x3 matrix
      {
         OMT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
         mat(0,0) = 1;
         mat(0,1) = 0;
         mat(0,2) = 2;
         mat(1,0) = 0;
         mat(1,1) = 3;
         mat(1,2) = 0;
         mat(2,0) = 4;
         mat(2,1) = 0;
         mat(2,2) = 5;

         transpose( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  5UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );
         checkNonZeros( mat,  2UL, 2UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 4 ||
             mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 ||
             mat(2,0) != 2 || mat(2,1) != 0 || mat(2,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 4 )\n( 0 3 0 )\n( 2 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Try to self-transpose a 5x3 matrix
      try {
         OMT mat( blaze::allocate<int>( 48UL ), 5UL, 3UL, 16UL, blaze::Deallocate() );

         transpose( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Self-transpose of a non-square matrix succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::logic_error& ) {}
   }

   {
      test_ = "Column-major self-transpose via trans()";

      // Self-transpose of a 3x3 matrix
      {
         OMT mat( blaze::allocate<int>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
         mat(0,0) = 1;
         mat(0,1) = 0;
         mat(0,2) = 2;
         mat(1,0) = 0;
         mat(1,1) = 3;
         mat(1,2) = 0;
         mat(2,0) = 4;
         mat(2,1) = 0;
         mat(2,2) = 5;

         mat = trans( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  5UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );
         checkNonZeros( mat,  2UL, 2UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 4 ||
             mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 ||
             mat(2,0) != 2 || mat(2,1) != 0 || mat(2,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 4 )\n( 0 3 0 )\n( 2 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Try to self-transpose a 5x3 matrix
      try {
         OMT mat( blaze::allocate<int>( 48UL ), 5UL, 3UL, 16UL, blaze::Deallocate() );

         mat = trans( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Self-transpose of a non-square matrix succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the CustomMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the CustomMatrix
// class template. Additionally, it performs a test of self-transpose via the \c ctrans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testCTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via ctranspose()";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::complex<int>  cplx;
      typedef blaze::CustomMatrix<cplx,aligned,padded,rowMajor>  AlignedPadded;

      // Self-transpose of a 3x3 matrix
      {
         AlignedPadded mat( blaze::allocate<cplx>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
         mat(0,0) = cplx(1,-1);
         mat(0,1) = cplx(0, 0);
         mat(0,2) = cplx(2,-2);
         mat(1,0) = cplx(0, 0);
         mat(1,1) = cplx(3,-3);
         mat(1,2) = cplx(0, 0);
         mat(2,0) = cplx(4,-4);
         mat(2,1) = cplx(0, 0);
         mat(2,2) = cplx(5,-5);

         ctranspose( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  5UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );
         checkNonZeros( mat,  2UL, 2UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(4,4) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(3,3) || mat(1,2) != cplx(0,0) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(5,5) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (4,4) )\n"
                                        "( (0,0) (3,3) (0,0) )\n"
                                        "( (2,2) (0,0) (5,5) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Try to self-transpose a 3x5 matrix
      try {
         AlignedPadded mat( blaze::allocate<cplx>( 48UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );

         ctranspose( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Self-transpose of a non-square matrix succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::logic_error& ) {}
   }

   {
      test_ = "Row-major self-transpose via ctrans()";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      typedef blaze::complex<int>  cplx;
      typedef blaze::CustomMatrix<cplx,aligned,padded,rowMajor>  AlignedPadded;

      // Self-transpose of a 3x3 matrix
      {
         AlignedPadded mat( blaze::allocate<cplx>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
         mat(0,0) = cplx(1,-1);
         mat(0,1) = cplx(0, 0);
         mat(0,2) = cplx(2,-2);
         mat(1,0) = cplx(0, 0);
         mat(1,1) = cplx(3,-3);
         mat(1,2) = cplx(0, 0);
         mat(2,0) = cplx(4,-4);
         mat(2,1) = cplx(0, 0);
         mat(2,2) = cplx(5,-5);

         mat = ctrans( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  5UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );
         checkNonZeros( mat,  2UL, 2UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(4,4) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(3,3) || mat(1,2) != cplx(0,0) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(5,5) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (4,4) )\n"
                                        "( (0,0) (3,3) (0,0) )\n"
                                        "( (2,2) (0,0) (5,5) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Try to self-transpose a 3x5 matrix
      try {
         AlignedPadded mat( blaze::allocate<cplx>( 48UL ), 3UL, 5UL, 16UL, blaze::Deallocate() );

         mat = ctrans( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Self-transpose of a non-square matrix succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::complex<int>  cplx;
      typedef blaze::CustomMatrix<cplx,aligned,padded,columnMajor>  AlignedPadded;

      // Self-transpose of a 3x3 matrix
      {
         AlignedPadded mat( blaze::allocate<cplx>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
         mat(0,0) = cplx(1,-1);
         mat(0,1) = cplx(0, 0);
         mat(0,2) = cplx(2,-2);
         mat(1,0) = cplx(0, 0);
         mat(1,1) = cplx(3,-3);
         mat(1,2) = cplx(0, 0);
         mat(2,0) = cplx(4,-4);
         mat(2,1) = cplx(0, 0);
         mat(2,2) = cplx(5,-5);

         ctranspose( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  5UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );
         checkNonZeros( mat,  2UL, 2UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(4,4) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(3,3) || mat(1,2) != cplx(0,0) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(5,5) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (4,4) )\n"
                                        "( (0,0) (3,3) (0,0) )\n"
                                        "( (2,2) (0,0) (5,5) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Try to self-transpose a 5x3 matrix
      try {
         AlignedPadded mat( blaze::allocate<cplx>( 48UL ), 5UL, 3UL, 16UL, blaze::Deallocate() );

         ctranspose( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Self-transpose of a non-square matrix succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::logic_error& ) {}
   }

   {
      test_ = "Column-major self-transpose via ctrans()";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      typedef blaze::complex<int>  cplx;
      typedef blaze::CustomMatrix<cplx,aligned,padded,columnMajor>  AlignedPadded;

      // Self-transpose of a 3x3 matrix
      {
         AlignedPadded mat( blaze::allocate<cplx>( 48UL ), 3UL, 3UL, 16UL, blaze::Deallocate() );
         mat(0,0) = cplx(1,-1);
         mat(0,1) = cplx(0, 0);
         mat(0,2) = cplx(2,-2);
         mat(1,0) = cplx(0, 0);
         mat(1,1) = cplx(3,-3);
         mat(1,2) = cplx(0, 0);
         mat(2,0) = cplx(4,-4);
         mat(2,1) = cplx(0, 0);
         mat(2,2) = cplx(5,-5);

         mat = ctrans( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
         checkNonZeros( mat,  5UL );
         checkNonZeros( mat,  0UL, 2UL );
         checkNonZeros( mat,  1UL, 1UL );
         checkNonZeros( mat,  2UL, 2UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(4,4) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(3,3) || mat(1,2) != cplx(0,0) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(5,5) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (4,4) )\n"
                                        "( (0,0) (3,3) (0,0) )\n"
                                        "( (2,2) (0,0) (5,5) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Try to self-transpose a 5x3 matrix
      try {
         AlignedPadded mat( blaze::allocate<cplx>( 48UL ), 5UL, 3UL, 16UL, blaze::Deallocate() );

         mat = ctrans( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Self-transpose of a non-square matrix succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the CustomMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the CustomMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix swap";

      MT mat1( blaze::allocate<int>( 32UL ), 2UL, 2UL, 16UL, blaze::Deallocate() );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 0;
      mat1(1,1) = 3;

      MT mat2( blaze::allocate<int>( 64UL ), 2UL, 2UL, 32UL, blaze::Deallocate() );
      mat2(0,0) = 4;
      mat2(0,1) = 3;
      mat2(1,0) = 2;
      mat2(1,1) = 1;

      swap( mat1, mat2 );

      checkRows    ( mat1,  2UL );
      checkColumns ( mat1,  2UL );
      checkCapacity( mat1, 64UL );
      checkNonZeros( mat1,  4UL );
      checkNonZeros( mat1,  0UL, 2UL );
      checkNonZeros( mat1,  1UL, 2UL );

      if( mat1(0,0) != 4 || mat1(0,1) != 3 || mat1(1,0) != 2 || mat1(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n( 4 3 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  2UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  3UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(1,0) != 0 || mat2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix swap";

      OMT mat1( blaze::allocate<int>( 32UL ), 2UL, 2UL, 16UL, blaze::Deallocate() );
      mat1(0,0) = 1;
      mat1(0,1) = 0;
      mat1(1,0) = 2;
      mat1(1,1) = 3;

      OMT mat2( blaze::allocate<int>( 64UL ), 2UL, 2UL, 32UL, blaze::Deallocate() );
      mat2(0,0) = 4;
      mat2(0,1) = 2;
      mat2(1,0) = 3;
      mat2(1,1) = 1;

      swap( mat1, mat2 );

      checkRows    ( mat1,  2UL );
      checkColumns ( mat1,  2UL );
      checkCapacity( mat1, 64UL );
      checkNonZeros( mat1,  4UL );
      checkNonZeros( mat1,  0UL, 2UL );
      checkNonZeros( mat1,  1UL, 2UL );

      if( mat1(0,0) != 4 || mat1(0,1) != 2 || mat1(1,0) != 3 || mat1(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n( 4 2 )\n( 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  2UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  3UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(1,0) != 2 || mat2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the CustomMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the CustomMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedPaddedTest::testIsDefault()
{
   using blaze::isDefault;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         MT mat;

         if( isDefault( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         MT mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
         reset( mat );

         if( isDefault( mat(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << mat(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         MT mat( blaze::allocate<int>( 64UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
         reset( mat );
         mat(0,1) = 1;

         if( isDefault( mat(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << mat(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         OMT mat;

         if( isDefault( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
         reset( mat );

         if( isDefault( mat(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << mat(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         OMT mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
         reset( mat );
         mat(1,0) = 1;

         if( isDefault( mat(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << mat(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************

} // namespace custommatrix

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
   std::cout << "   Running aligned/padded CustomMatrix class test..." << std::endl;

   try
   {
      RUN_CUSTOMMATRIX_ALIGNED_PADDED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aligned/padded CustomMatrix class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
