//=================================================================================================
/*!
//  \file src/mathtest/unilowermatrix/SparseTest.cpp
//  \brief Source file for the UniLowerMatrix sparse test
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/SparseColumn.h>
#include <blaze/math/SparseRow.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/unilowermatrix/SparseTest.h>


namespace blazetest {

namespace mathtest {

namespace unilowermatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the UniLowerMatrix sparse test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseTest::SparseTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testFunctionCall();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testSet();
   testInsert();
   testAppend();
   testErase();
   testResize();
   testReserve();
   testTrim();
   testSwap();
   testFind();
   testLowerBound();
   testUpperBound();
   testIsDefault();
   testSubmatrix();
   testRow();
   testColumn();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the UniLowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (CompressedMatrix)
   {
      test_ = "Row-major UniLowerMatrix default constructor (CompressedMatrix)";

      const LT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   // Size constructor (CompressedMatrix)
   {
      test_ = "Row-major UniLowerMatrix size constructor (CompressedMatrix)";

      const LT lower( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkNonZeros( lower, 2UL );
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Row-major UniLowerMatrix copy constructor (0x0)";

      const LT lower1;
      const LT lower2( lower1 );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Row-major UniLowerMatrix copy constructor (3x3)";

      LT lower1( 3UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      const LT lower2( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Row-major UniLowerMatrix move constructor (0x0)";

      LT lower1;
      LT lower2( std::move( lower1 ) );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Row-major UniLowerMatrix move constructor (3x3)";

      LT lower1( 3UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      LT lower2( std::move( lower1 ) );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Row-major UniLowerMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      const LT lower( mat );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Conversion constructor (unilower)
   {
      test_ = "Row-major UniLowerMatrix conversion constructor (unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      const LT lower( mat );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-unilower)
   {
      test_ = "Row-major UniLowerMatrix conversion constructor (non-unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         const LT lower( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-unilower UniLowerMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (UniLowerMatrix)
   {
      test_ = "Row-major UniLowerMatrix conversion constructor (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      const LT lower2( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (CompressedMatrix)
   {
      test_ = "Column-major UniLowerMatrix default constructor (CompressedMatrix)";

      const OLT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   // Size constructor (CompressedMatrix)
   {
      test_ = "Column-major UniLowerMatrix size constructor (CompressedMatrix)";

      const OLT lower( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkNonZeros( lower, 2UL );
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Column-major UniLowerMatrix copy constructor (0x0)";

      const OLT lower1;
      const OLT lower2( lower1 );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Column-major UniLowerMatrix copy constructor (3x3)";

      OLT lower1( 3UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      const OLT lower2( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Column-major UniLowerMatrix move constructor (0x0)";

      OLT lower1;
      OLT lower2( std::move( lower1 ) );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Column-major UniLowerMatrix move constructor (3x3)";

      OLT lower1( 3UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      OLT lower2( std::move( lower1 ) );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Column-major UniLowerMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      const OLT lower( mat );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Conversion constructor (unilower)
   {
      test_ = "Column-major UniLowerMatrix conversion constructor (unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      const OLT lower( mat );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-unilower)
   {
      test_ = "Column-major UniLowerMatrix conversion constructor (non-unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         const OLT lower( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-unilower UniLowerMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (UniLowerMatrix)
   {
      test_ = "Column-major UniLowerMatrix conversion constructor (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      const OLT lower2( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the UniLowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAssignment()
{
   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Row-major UniLowerMatrix copy assignment (0x0)";

      LT lower1, lower2;

      lower2 = lower1;

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Row-major UniLowerMatrix copy assignment (3x3)";

      LT lower1( 3UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Row-major UniLowerMatrix move assignment (0x0)";

      LT lower1, lower2;

      lower2 = std::move( lower1 );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Row-major UniLowerMatrix move assignment (3x3)";

      LT lower1( 3UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      LT lower2;
      lower2 = std::move( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major UniLowerMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      LT lower;
      lower = mat;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Row-major/row-major dense matrix assignment (unilower)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix assignment (unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      LT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (unilower)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix assignment (unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      LT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (non-unilower)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix assignment (non-unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      try {
         LT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-unilower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (non-unilower)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix assignment (non-unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      try {
         LT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-unilower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (UniLowerMatrix)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (UniLowerMatrix)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major UniLowerMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      LT lower;
      lower = mat;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (unilower)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix assignment (unilower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;
      mat.insert( 1UL, 2UL, 0 );

      LT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (unilower)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix assignment (unilower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;
      mat.insert( 1UL, 2UL, 0 );

      LT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (non-unilower)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix assignment (non-unilower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      try {
         LT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-unilower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (non-unilower)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix assignment (non-unilower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      try {
         LT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-unilower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (UniLowerMatrix)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::CompressedMatrix<unsigned int,blaze::rowMajor> > lower1( 3UL, 5UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (UniLowerMatrix)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > lower1( 3UL, 5UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Column-major UniLowerMatrix copy assignment (0x0)";

      OLT lower1, lower2;

      lower2 = lower1;

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Column-major UniLowerMatrix copy assignment (3x3)";

      OLT lower1( 3UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Column-major UniLowerMatrix move assignment (0x0)";

      OLT lower1, lower2;

      lower2 = std::move( lower1 );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Column-major UniLowerMatrix move assignment (3x3)";

      OLT lower1( 3UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      OLT lower2;
      lower2 = std::move( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major UniLowerMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      OLT lower;
      lower = mat;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Column-major/row-major dense matrix assignment (unilower)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix assignment (unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      OLT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (unilower)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix assignment (unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      OLT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (non-unilower)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix assignment (non-unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      try {
         OLT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-unilower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (non-unilower)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix assignment (non-unilower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      try {
         OLT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-unilower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (UniLowerMatrix)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (UniLowerMatrix)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major UniLowerMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat;

      OLT lower;
      lower = mat;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (unilower)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix assignment (unilower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (unilower)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix assignment (unilower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (non-unilower)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix assignment (non-unilower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      try {
         OLT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-unilower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (non-unilower)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix assignment (non-unilower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  1;
      mat(2,0) =  7;
      mat(2,2) =  1;

      try {
         OLT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-unilower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (UniLowerMatrix)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > lower1( 3UL, 5UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (UniLowerMatrix)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::CompressedMatrix<unsigned int,blaze::columnMajor> > lower1( 3UL, 5UL );
      lower1(1,0) = -4;
      lower1(2,0) =  7;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   // Row-major/row-major dense matrix addition assignment (strictly lower)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix addition assignment (strictly lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) =  2;
      mat(2,0) = -7;
      mat(2,1) =  5;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (strictly lower)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix addition assignment (strictly lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) =  2;
      mat(2,0) = -7;
      mat(2,1) =  5;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix addition assignment (non-lower)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix addition assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,2) = 6;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (non-lower)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix addition assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,2) = 6;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix addition assignment (strictly lower)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix addition assignment (strictly lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,0) =  2;
      mat(2,0) = -7;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (strictly lower)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix addition assignment (strictly lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,0) =  2;
      mat(2,0) = -7;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix addition assignment (non-lower)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix addition assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,2) = 6;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix addition assignment (non-lower)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix addition assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,2) = 6;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   // Column-major/row-major dense matrix addition assignment (strictly lower)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix addition assignment (strictly lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) =  2;
      mat(2,0) = -7;
      mat(2,1) =  5;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (strictly lower)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix addition assignment (strictly lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) =  2;
      mat(2,0) = -7;
      mat(2,1) =  5;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix addition assignment (non-lower)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix addition assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,2) = 6;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (non-lower)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix addition assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,2) = 6;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix addition assignment (strictly lower)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix addition assignment (strictly lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,0) =  2;
      mat(2,0) = -7;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (strictly lower)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix addition assignment (strictly lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,0) =  2;
      mat(2,0) = -7;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix addition assignment (non-lower)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix addition assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,2) = 6;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix addition assignment (non-lower)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix addition assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,2) = 6;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major dense matrix subtraction assignment (strictly lower)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix subtraction assignment (strictly lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) = -2;
      mat(2,0) =  7;
      mat(2,1) =  5;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  0 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (strictly lower)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix subtraction assignment (strictly lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) = -2;
      mat(2,0) =  7;
      mat(2,1) =  5;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  1 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix subtraction assignment (non-lower)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix subtraction assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,2) = 6;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (non-lower)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix subtraction assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,2) = 6;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix subtraction assignment (strictly lower)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix subtraction assignment (strictly lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,0) = -2;
      mat(2,0) =  7;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  0 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (strictly lower)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix subtraction assignment (strictly lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,0) = -2;
      mat(2,0) =  7;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  0 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix subtraction assignment (non-lower)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix subtraction assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,2) = 6;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix subtraction assignment (non-lower)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix subtraction assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,2) = 6;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major dense matrix subtraction assignment (strictly lower)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix subtraction assignment (strictly lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) = -2;
      mat(2,0) =  7;
      mat(2,1) =  5;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  0 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (strictly lower)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix subtraction assignment (strictly lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) = -2;
      mat(2,0) =  7;
      mat(2,1) =  5;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  0 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix subtraction assignment (non-lower)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix subtraction assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,2) = 6;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (non-lower)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix subtraction assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,2) = 6;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix subtraction assignment (strictly lower)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix subtraction assignment (strictly lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,0) = -2;
      mat(2,0) =  7;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  0 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (strictly lower)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix subtraction assignment (strictly lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,0) = -2;
      mat(2,0) =  7;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  0 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  0 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix subtraction assignment (non-lower)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix subtraction assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,2) = 6;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix subtraction assignment (non-lower)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix subtraction assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,2) = 6;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major dense matrix multiplication assignment (unilower)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix multiplication assignment (unilower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) =  1;
      mat(1,1) =  1;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  5 || lower(2,1) != 3 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (unilower)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix multiplication assignment (unilower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) =  1;
      mat(1,1) =  1;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  5 || lower(2,1) != 3 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix multiplication assignment (non-unilower)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix multiplication assignment (non-unilower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) =  1;
      mat(1,1) =  4;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-unilower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (non-unilower)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix multiplication assignment (non-unilower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) =  1;
      mat(1,1) =  4;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-unilower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (UniLowerMatrix)
   {
      test_ = "Row-major/row-major UniLowerMatrix dense matrix multiplication assignment (UniLowerMatrix)";

      LT lower1( 3UL );
      lower1(2,0) = -2;
      lower1(2,1) =  3;

      LT lower2( 3UL );
      lower2(1,0) = -4;
      lower2(2,0) =  7;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  5 || lower2(2,1) != 3 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (UniLowerMatrix)
   {
      test_ = "Row-major/column-major UniLowerMatrix dense matrix multiplication assignment (UniLowerMatrix)";

      OLT lower1( 3UL );
      lower1(2,0) = -2;
      lower1(2,1) =  3;

      LT lower2( 3UL );
      lower2(1,0) = -4;
      lower2(2,0) =  7;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  5 || lower2(2,1) != 3 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix multiplication assignment (unilower)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix multiplication assignment (unilower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,1) =  1;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  5 || lower(2,1) != 3 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (unilower)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix multiplication assignment (unilower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,1) =  1;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  5 || lower(2,1) != 3 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix multiplication assignment (non-unilower)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix multiplication assignment (non-unilower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) =  1;
      mat(1,1) =  4;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-unilower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (non-unilower)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix multiplication assignment (non-unilower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) =  1;
      mat(1,1) =  4;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-unilower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (UniLowerMatrix)
   {
      test_ = "Row-major/row-major UniLowerMatrix sparse matrix multiplication assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > lower1( 3UL, 5UL );
      lower1(2,0) = -2;
      lower1(2,1) =  3;

      LT lower2( 3UL );
      lower2(1,0) = -4;
      lower2(2,0) =  7;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  5 || lower2(2,1) != 3 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (UniLowerMatrix)
   {
      test_ = "Row-major/column-major UniLowerMatrix sparse matrix multiplication assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > lower1( 3UL, 5UL );
      lower1(2,0) = -2;
      lower1(2,1) =  3;

      LT lower2( 3UL );
      lower2(1,0) = -4;
      lower2(2,0) =  7;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  5 || lower2(2,1) != 3 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major dense matrix multiplication assignment (unilower)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix multiplication assignment (unilower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) =  1;
      mat(1,1) =  1;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  5 || lower(2,1) != 3 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (unilower)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix multiplication assignment (unilower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) =  1;
      mat(1,1) =  1;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  5 || lower(2,1) != 3 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix multiplication assignment (non-unilower)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix multiplication assignment (non-unilower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) =  1;
      mat(1,1) =  4;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-unilower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (non-unilower)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix multiplication assignment (non-unilower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) =  1;
      mat(1,1) =  4;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-unilower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (UniLowerMatrix)
   {
      test_ = "Column-major/row-major UniLowerMatrix dense matrix multiplication assignment (UniLowerMatrix)";

      LT lower1( 3UL );
      lower1(2,0) = -2;
      lower1(2,1) =  3;

      OLT lower2( 3UL );
      lower2(1,0) = -4;
      lower2(2,0) =  7;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  5 || lower2(2,1) != 3 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (UniLowerMatrix)
   {
      test_ = "Column-major/column-major UniLowerMatrix dense matrix multiplication assignment (UniLowerMatrix)";

      OLT lower1( 3UL );
      lower1(2,0) = -2;
      lower1(2,1) =  3;

      OLT lower2( 3UL );
      lower2(1,0) = -4;
      lower2(2,0) =  7;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  5 || lower2(2,1) != 3 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix multiplication assignment (unilower)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix multiplication assignment (unilower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,1) =  1;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  5 || lower(2,1) != 3 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (unilower)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix multiplication assignment (unilower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,1) =  1;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  5 || lower(2,1) != 3 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix multiplication assignment (non-unilower)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix multiplication assignment (non-unilower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) =  1;
      mat(1,1) =  4;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-unilower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (non-unilower)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix multiplication assignment (non-unilower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) =  1;
      mat(1,1) =  4;
      mat(2,0) = -2;
      mat(2,1) =  3;
      mat(2,2) =  1;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-unilower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (UniLowerMatrix)
   {
      test_ = "Column-major/row-major UniLowerMatrix sparse matrix multiplication assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > lower1( 3UL, 5UL );
      lower1(2,0) = -2;
      lower1(2,1) =  3;

      OLT lower2( 3UL );
      lower2(1,0) = -4;
      lower2(2,0) =  7;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  5 || lower2(2,1) != 3 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (UniLowerMatrix)
   {
      test_ = "Column-major/column-major UniLowerMatrix sparse matrix multiplication assignment (UniLowerMatrix)";

      blaze::UniLowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > lower1( 3UL, 5UL );
      lower1(2,0) = -2;
      lower1(2,1) =  3;

      OLT lower2( 3UL );
      lower2(1,0) = -4;
      lower2(2,0) =  7;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 1 || lower2(1,2) != 0 ||
          lower2(2,0) !=  5 || lower2(2,1) != 3 || lower2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  5 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the UniLowerMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::operator()";

      // Good cases
      {
         LT lower( 3UL );

         // Writing the element (2,1)
         lower(2,1) = 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the element (1,0)
         lower(1,0) = lower(2,1);

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the element (2,0)
         lower(2,0) += 3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 3 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the element (1,0)
         lower(1,0) -= 4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  3 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the element (2,1)
         lower(2,1) *= -3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -6 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  3 -6  1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the element (2,1)
         lower(2,1) /= 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -3 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  3 -3  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         LT lower( 3UL );

         // Trying to write the diagonal element (1,1)
         try {
            lower(1,1) = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (1,2)
         try {
            lower(1,2) = 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the diagonal element (2,2)
         try {
            lower(2,2) = lower(1,1);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (0,1)
         try {
            lower(0,1) = lower(2,1);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the diagonal element (1,1)
         try {
            lower(1,1) += 6;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the upper element (0,2)
         try {
            lower(0,2) += 3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the diagonal element (1,1)
         try {
            lower(1,1) -= 8;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the upper element (0,1)
         try {
            lower(0,1) -= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the diagonal element (1,1)
         try {
            lower(1,1) *= -6;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the upper element (1,2)
         try {
            lower(1,2) *= -3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the diagonal element (1,1)
         try {
            lower(1,1) /= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the upper element (1,2)
         try {
            lower(1,2) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::operator()";

      // Good cases
      {
         OLT lower( 3UL );

         // Writing the lower element (2,1)
         lower(2,1) = 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (1,0)
         lower(1,0) = lower(2,1);

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the lower element (2,0)
         lower(2,0) += 3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 3 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the lower element (1,0)
         lower(1,0) -= 4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  3 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the lower element (2,1)
         lower(2,1) *= -3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -6 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  3 -6  1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the lower element (2,1)
         lower(2,1) /= 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -3 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  3 -3  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         OLT lower( 3UL );

         // Trying to write the diagonal element (1,1)
         try {
            lower(1,1) = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (1,2)
         try {
            lower(1,2) = 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the diagonal element (2,2)
         try {
            lower(2,2) = lower(1,1);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (0,1)
         try {
            lower(0,1) = lower(2,1);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the diagonal element (1,1)
         try {
            lower(1,1) += 6;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the upper element (0,2)
         try {
            lower(0,2) += 3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the diagonal element (1,1)
         try {
            lower(1,1) -= 8;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the upper element (0,1)
         try {
            lower(0,1) -= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the diagonal element (1,1)
         try {
            lower(1,1) *= -6;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the upper element (1,2)
         try {
            lower(1,2) *= -3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the diagonal element (1,1)
         try {
            lower(1,1) /= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the upper element (1,2)
         try {
            lower(1,2) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      typedef LT::Iterator       Iterator;
      typedef LT::ConstIterator  ConstIterator;

      LT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 3;

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

         ConstIterator it( begin( lower, 1UL ) );

         if( it == end( lower, 1UL ) || it->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator
      {
         test_ = "Row-major Iterator subtraction";

         const size_t number( end( lower, 0UL ) - begin( lower, 0UL ) );

         if( number != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st row via ConstIterator
      {
         test_ = "Row-major ConstIterator subtraction";

         const size_t number( cend( lower, 1UL ) - cbegin( lower, 1UL ) );

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
         test_ = "Row-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( lower, 2UL ) );
         ConstIterator end( cend( lower, 2UL ) );

         if( it == end || it->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != 1 ) {
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

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Row-major assignment to lower elements via Iterator";

         Iterator it1 = begin( lower, 1UL );
         Iterator it2 = begin( lower, 2UL );
         *it1 = 5;
         it2->value() = 7;

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 5 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 5 1 0 )\n( 7 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to diagonal elements via Iterator
      {
         test_ = "Row-major assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Row-major addition assignment to lower elements via Iterator";

         Iterator it1 = begin( lower, 1UL );
         Iterator it2 = begin( lower, 2UL );
         *it1 += 2;
         it2->value() += -2;

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 7 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 5 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 7 1 0 )\n( 5 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to diagonal elements via Iterator
      {
         test_ = "Row-major addition assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Row-major subtraction assignment to lower elements via Iterator";

         Iterator it1 = begin( lower, 1UL );
         Iterator it2 = begin( lower, 2UL );
         *it1 -= 2;
         it2->value() -= -2;

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 5 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 5 1 0 )\n( 7 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to diagonal elements via Iterator
      {
         test_ = "Row-major subtraction assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Row-major multiplication assignment to lower elements via Iterator";

         Iterator it1 = begin( lower, 1UL );
         Iterator it2 = begin( lower, 2UL );
         *it1 *= 2;
         it2->value() *= -2;

         if( lower(0,0) !=   1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) !=  10 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != -14 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(   1 0 0 )\n(  10 1 0 )\n( -14 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to diagonal elements via Iterator
      {
         test_ = "Row-major multiplication assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Row-major division assignment to lower elements via Iterator";

         Iterator it1 = begin( lower, 1UL );
         Iterator it2 = begin( lower, 2UL );
         *it1 /= 2;
         it2->value() /= -2;

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 5 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 5 1 0 )\n( 7 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to diagonal elements via Iterator
      {
         test_ = "Row-major division assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      typedef OLT::Iterator       Iterator;
      typedef OLT::ConstIterator  ConstIterator;

      OLT lower( 3UL );
      lower(2,0) = 3;
      lower(2,1) = 2;

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

         ConstIterator it( begin( lower, 1UL ) );

         if( it == end( lower, 1UL ) || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator
      {
         test_ = "Column-major Iterator subtraction";

         const size_t number( end( lower, 0UL ) - begin( lower, 0UL ) );

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

      // Counting the number of elements in 1st column via ConstIterator
      {
         test_ = "Column-major ConstIterator subtraction";

         const size_t number( cend( lower, 1UL ) - cbegin( lower, 1UL ) );

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
         test_ = "Column-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( lower, 0UL ) );
         ConstIterator end( cend( lower, 0UL ) );

         if( it == end || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != 3 ) {
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

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Column-major assignment to lower elements via Iterator";

         Iterator it1 = lower.find( 2UL, 0UL );
         Iterator it2 = lower.find( 2UL, 1UL );
         *it1 = 5;
         it2->value() = 7;

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 5 || lower(2,1) != 7 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 5 7 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to diagonal elements via Iterator
      {
         test_ = "Column-major assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Column-major addition assignment to lower elements via Iterator";

         Iterator it1 = lower.find( 2UL, 0UL );
         Iterator it2 = lower.find( 2UL, 1UL );
         *it1 += 2;
         it2->value() += -2;

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 7 || lower(2,1) != 5 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 7 5 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to diagonal elements via Iterator
      {
         test_ = "Column-major addition assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Column-major subtraction assignment to lower elements via Iterator";

         Iterator it1 = lower.find( 2UL, 0UL );
         Iterator it2 = lower.find( 2UL, 1UL );
         *it1 -= 2;
         it2->value() -= -2;

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 5 || lower(2,1) != 7 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 5 7 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to diagonal elements via Iterator
      {
         test_ = "Column-major subtraction assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Column-major multiplication assignment to lower elements via Iterator";

         Iterator it1 = lower.find( 2UL, 0UL );
         Iterator it2 = lower.find( 2UL, 1UL );
         *it1 *= 2;
         it2->value() *= -2;

         if( lower(0,0) !=  1 || lower(0,1) !=   0 || lower(0,2) != 0 ||
             lower(1,0) !=  0 || lower(1,1) !=   1 || lower(1,2) != 0 ||
             lower(2,0) != 10 || lower(2,1) != -14 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1   0  0 )\n(  0   1  0 )\n( 10 -14  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to diagonal elements via Iterator
      {
         test_ = "Column-major multiplication assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Column-major division assignment to lower elements via Iterator";

         Iterator it1 = lower.find( 2UL, 0UL );
         Iterator it2 = lower.find( 2UL, 1UL );
         *it1 /= 2;
         it2->value() /= -2;

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 5 || lower(2,1) != 7 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 5 7 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to diagonal elements via Iterator
      {
         test_ = "Column-major division assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         try {
            const Iterator it = begin( lower, 0UL );
            it->value() /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::nonZeros()";

      // Default matrix
      {
         LT lower( 3UL );

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         LT lower( 3UL );
         lower(1,0) =  2;
         lower(2,0) = -4;
         lower(2,1) = -5;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != -4 || lower(2,1) != -5 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n(  2  1  0 )\n( -4 -5  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::nonZeros()";

      // Default matrix
      {
         OLT lower( 3UL );

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         OLT lower( 3UL );
         lower(1,0) =  2;
         lower(2,0) = -4;
         lower(2,1) = -5;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != -4 || lower(2,1) != -5 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n(  2  1  0 )\n( -4 -5  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::reset()";

      // Initialization check
      LT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a diagonal element
      reset( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( lower, 1UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::reset()";

      // Initialization check
      OLT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a diagonal element
      reset( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting column 1
      reset( lower, 1UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::clear()";

      // Initialization check
      LT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( lower );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::clear()";

      // Initialization check
      OLT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( lower );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::set()";

      typedef LT::Iterator  Iterator;

      // Initialization check
      LT lower( 4UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      // Setting a non-zero element
      {
         Iterator pos = lower.set( 2UL, 1UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 2 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 2 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = lower.set( 2UL, 0UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 3UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 3 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 3 2 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = lower.set( 2UL, 1UL, 4 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 3UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 4 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 3 4 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::set()";

      typedef OLT::Iterator  Iterator;

      // Initialization check
      OLT lower( 4UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      // Setting a non-zero element
      {
         Iterator pos = lower.set( 2UL, 1UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 2 || pos->index() != 2UL ) {
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

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 2 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = lower.set( 3UL, 1UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 3 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,1) != 3 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 2 1 0 )\n( 0 3 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = lower.set( 2UL, 1UL, 4 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 4 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,1) != 4 || lower(2,2) != 1 ||
             lower(3,1) != 3 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 4 1 0 )\n( 0 3 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::insert()";

      typedef LT::Iterator  Iterator;

      // Initialization check
      LT lower( 4UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      // Inserting a non-zero element
      {
         Iterator pos = lower.insert( 2UL, 1UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 2 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 2 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = lower.insert( 2UL, 0UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 3UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 3 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 3 2 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         lower.insert( 2UL, 1UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 3 2 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::set()";

      typedef OLT::Iterator  Iterator;

      // Initialization check
      OLT lower( 4UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      // Inserting a non-zero element
      {
         Iterator pos = lower.set( 3UL, 1UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( pos->value() != 2 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,2) != 1 ||
             lower(3,1) != 2 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 2 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = lower.set( 2UL, 1UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

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

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,1) != 3 || lower(2,2) != 1 ||
             lower(3,1) != 2 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 3 1 0 )\n( 0 2 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         lower.insert( 2UL, 1UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 3 1 0 )\n( 0 2 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::append()";

      // Initialization check
      LT lower( 4UL, 5UL );
      lower.reserve( 2UL, 2UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      // Trying to append an element
      try {
         lower.append( 2UL, 3UL, 2 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Appending an upper element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::append()";

      // Appending with pre-allocation in each column
      {
         // Initialization check
         OLT lower( 4UL, 9UL );
         lower.reserve( 0UL, 3UL );
         lower.reserve( 1UL, 3UL );
         lower.reserve( 2UL, 2UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         // Appending one non-zero element
         lower.append( 2UL, 1UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 2 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 1UL, 0UL, 3 );
         lower.append( 3UL, 2UL, 4 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 7UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 3 || lower(1,1) != 1 ||
             lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,2) != 4 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 3 1 0 0 )\n( 0 2 1 0 )\n( 0 0 4 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 3UL, 0UL, 5 );
         lower.append( 3UL, 1UL, 6 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 9UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 3 || lower(1,1) != 1 ||
             lower(2,1) != 2 || lower(2,2) != 1 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,2) != 4 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 3 1 0 0 )\n( 0 2 1 0 )\n( 5 6 4 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with column finalization
      {
         // Initialization check
         OLT lower( 4UL, 8UL );
         lower.reserve( 0UL, 2UL );
         lower.reserve( 1UL, 3UL );
         lower.reserve( 2UL, 2UL );

         // Appending one non-zero element
         lower.append( 1UL, 0UL, 2 );
         lower.finalize( 0UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 2UL, 1UL, 3 );
         lower.append( 3UL, 1UL, 4 );
         lower.finalize( 1UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 7UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,1) != 3 || lower(2,2) != 1 ||
             lower(3,1) != 4 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 0 3 1 0 )\n( 0 4 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending one more non-zero element
         lower.append( 3UL, 2UL, 5 );
         lower.finalize( 2UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 8UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,1) != 3 || lower(2,2) != 1 ||
             lower(3,1) != 4 || lower(3,2) != 5 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 0 3 1 0 )\n( 0 4 5 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::erase( size_t, size_t )";

      // Initialization check
      LT lower( 4UL, 9UL );
      lower(1,0) = 2;
      lower(2,0) = 3;
      lower(2,1) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 4 1 0 )\n( 5 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      lower.erase( 2UL, 1UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,2) != 1 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 5 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (3,0)
      lower.erase( 3UL, 0UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 7UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 2UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,2) != 1 ||
          lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      lower.erase( 3UL, 2UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 7UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 2UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,2) != 1 ||
          lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a diagonal element
      try {
         lower.erase( 0UL, 0UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::erase( size_t, Iterator )";

      typedef LT::Iterator  Iterator;

      // Initialization check
      LT lower( 4UL, 9UL );
      lower(1,0) = 2;
      lower(2,0) = 3;
      lower(2,1) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 4 1 0 )\n( 5 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      {
         Iterator pos = lower.erase( 2UL, lower.find( 2UL, 1UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 8UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 3UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 5 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 1 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (3,0)
      {
         Iterator pos = lower.erase( 3UL, lower.find( 3UL, 0UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 2UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 6 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 6\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = lower.erase( 3UL, lower.find( 3UL, 2UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 2UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != lower.end( 3UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a diagonal element
      try {
         lower.erase( 0UL, lower.find( 0UL, 0UL ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::erase( size_t, Iterator, Iterator )";

      typedef LT::Iterator  Iterator;

      // Initialization check
      LT lower( 4UL, 9UL );
      lower(1,0) = 2;
      lower(2,0) = 3;
      lower(2,1) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 4 1 0 )\n( 5 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from (2,1) to (2,2)
      {
         Iterator pos = lower.erase( 2UL, lower.find( 2UL, 1UL ), lower.find( 2UL, 2UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 8UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 3UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 5 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 1 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from the beginning of row 3 to (3,3)
      {
         Iterator pos = lower.erase( 3UL, lower.begin( 3UL ), lower.find( 3UL, 3UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a multi-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 1 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = lower.erase( 3UL, lower.find( 3UL, 3UL ), lower.find( 3UL, 3UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 1 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a range including a diagonal element
      try {
         lower.erase( 1UL, lower.begin( 1UL ), lower.end( 1UL ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a range including a diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::erase( size_t, size_t )";

      // Initialization check
      OLT lower( 4UL, 9UL );
      lower(1,0) = 2;
      lower(2,0) = 3;
      lower(2,1) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 4UL );
      checkNonZeros( lower, 1UL, 3UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 4 1 0 )\n( 5 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      lower.erase( 2UL, 1UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 4UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,2) != 1 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 5 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (3,0)
      lower.erase( 3UL, 0UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,2) != 1 ||
          lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      lower.erase( 3UL, 2UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,2) != 1 ||
          lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a diagonal element
      try {
         lower.erase( 3UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::erase( size_t, Iterator )";

      typedef OLT::Iterator  Iterator;

      // Initialization check
      OLT lower( 4UL, 9UL );
      lower(1,0) = 2;
      lower(2,0) = 3;
      lower(2,1) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 4UL );
      checkNonZeros( lower, 1UL, 3UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 4 1 0 )\n( 5 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      {
         Iterator pos = lower.erase( 1UL, lower.find( 2UL, 1UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 8UL );
         checkNonZeros( lower, 0UL, 4UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 5 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 6 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 6\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (3,0)
      {
         Iterator pos = lower.erase( 0UL, lower.find( 3UL, 0UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != lower.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = lower.erase( 2UL, lower.find( 3UL, 2UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a diagonal element
      try {
         lower.erase( 3UL, lower.find( 3UL, 3UL ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::erase( size_t, Iterator, Iterator )";

      typedef OLT::Iterator  Iterator;

      // Initialization check
      OLT lower( 4UL, 9UL );
      lower(1,0) = 2;
      lower(2,0) = 3;
      lower(2,1) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 4UL );
      checkNonZeros( lower, 1UL, 3UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 4 1 0 )\n( 5 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from (2,1) to (3,1)
      {
         Iterator pos = lower.erase( 1UL, lower.find( 2UL, 1UL ), lower.find( 3UL, 1UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 8UL );
         checkNonZeros( lower, 0UL, 4UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,0) != 2 || lower(1,1) != 1 ||
             lower(2,0) != 3 || lower(2,2) != 1 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 3 0 1 0 )\n( 5 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 6 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 6\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (1,0) to the column end
      {
         Iterator pos = lower.erase( 0UL, lower.find( 1UL, 0UL ), lower.end( 0UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,2) != 1 ||
             lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a multi-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != lower.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = lower.erase( 3UL, lower.find( 3UL, 3UL ), lower.find( 3UL, 3UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 1 ||
             lower(2,2) != 1 ||
             lower(3,1) != 6 || lower(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 6 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 1 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a range including a diagonal element
      try {
         lower.erase( 2UL, lower.begin( 2UL ), lower.end( 2UL ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a range including a diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 6 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::resize()";

      // Initialization check
      LT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Resizing to 2x2
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 2UL );
      checkNonZeros( lower, 2UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );

      if( lower(0,0) != 1 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      lower(1,0) = 2;
      lower.resize( 4UL, true );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,2) != 1 ||
          lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      lower(2,1) = 4;
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      lower.resize( 0UL );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::resize()";

      // Initialization check
      OLT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Resizing to 2x2
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 2UL );
      checkNonZeros( lower, 2UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( x 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      lower(1,0) = 2;
      lower.resize( 4UL, true );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ||
          lower(2,2) != 1 ||
          lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      lower(2,1) = 4;
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      lower.resize( 0UL );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::reserve()";

      // Initialization check
      LT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Increasing the capacity of the matrix
      lower.reserve( 10UL );

      checkRows    ( lower,  0UL );
      checkColumns ( lower,  0UL );
      checkCapacity( lower, 10UL );
      checkNonZeros( lower,  0UL );

      // Further increasing the capacity of the matrix
      lower.reserve( 20UL );

      checkRows    ( lower,  0UL );
      checkColumns ( lower,  0UL );
      checkCapacity( lower, 20UL );
      checkNonZeros( lower,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::reserve()";

      // Initialization check
      OLT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Increasing the capacity of the matrix
      lower.reserve( 10UL );

      checkRows    ( lower,  0UL );
      checkColumns ( lower,  0UL );
      checkCapacity( lower, 10UL );
      checkNonZeros( lower,  0UL );

      // Further increasing the capacity of the matrix
      lower.reserve( 20UL );

      checkRows    ( lower,  0UL );
      checkColumns ( lower,  0UL );
      checkCapacity( lower, 20UL );
      checkNonZeros( lower,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trim() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::trim()";

      // Initialization check
      LT lower( 3UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      // Increasing the row capacity of the matrix
      lower.reserve( 0UL, 10UL );
      lower.reserve( 1UL, 15UL );
      lower.reserve( 2UL, 20UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 10UL );
      checkCapacity( lower,  1UL, 15UL );
      checkCapacity( lower,  2UL, 20UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );

      // Trimming the matrix
      lower.trim();

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 1UL );
      checkCapacity( lower,  1UL, 1UL );
      checkCapacity( lower,  2UL, 1UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );
   }

   {
      test_ = "Row-major UniLowerMatrix::trim( size_t )";

      // Initialization check
      LT lower( 3UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      // Increasing the row capacity of the matrix
      lower.reserve( 0UL, 10UL );
      lower.reserve( 1UL, 15UL );
      lower.reserve( 2UL, 20UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 10UL );
      checkCapacity( lower,  1UL, 15UL );
      checkCapacity( lower,  2UL, 20UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );

      // Trimming the 0th row
      lower.trim( 0UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL,  1UL );
      checkCapacity( lower,  1UL, 24UL );
      checkCapacity( lower,  2UL, 20UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );

      // Trimming the 1st row
      lower.trim( 1UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL,  1UL );
      checkCapacity( lower,  1UL,  1UL );
      checkCapacity( lower,  2UL, 43UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );

      // Trimming the 2nd row
      lower.trim( 2UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 1UL );
      checkCapacity( lower,  1UL, 1UL );
      checkCapacity( lower,  2UL, 1UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::trim()";

      // Initialization check
      OLT lower( 3UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      // Increasing the column capacity of the matrix
      lower.reserve( 0UL, 10UL );
      lower.reserve( 1UL, 15UL );
      lower.reserve( 2UL, 20UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 10UL );
      checkCapacity( lower,  1UL, 15UL );
      checkCapacity( lower,  2UL, 20UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );

      // Trimming the matrix
      lower.trim();

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 1UL );
      checkCapacity( lower,  1UL, 1UL );
      checkCapacity( lower,  2UL, 1UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );
   }

   {
      test_ = "Column-major UniLowerMatrix::trim( size_t )";

      // Initialization check
      OLT lower( 3UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      // Increasing the column capacity of the matrix
      lower.reserve( 0UL, 10UL );
      lower.reserve( 1UL, 15UL );
      lower.reserve( 2UL, 20UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 10UL );
      checkCapacity( lower,  1UL, 15UL );
      checkCapacity( lower,  2UL, 20UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );

      // Trimming the 0th column
      lower.trim( 0UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL,  1UL );
      checkCapacity( lower,  1UL, 24UL );
      checkCapacity( lower,  2UL, 20UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );

      // Trimming the 1st column
      lower.trim( 1UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL,  1UL );
      checkCapacity( lower,  1UL,  1UL );
      checkCapacity( lower,  2UL, 43UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );

      // Trimming the 2nd column
      lower.trim( 2UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 1UL );
      checkCapacity( lower,  1UL, 1UL );
      checkCapacity( lower,  2UL, 1UL );
      checkNonZeros( lower,  3UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 1UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the UniLowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix swap";

      LT lower1( 2UL );
      lower1(1,0) = 2;

      LT lower2( 3UL );
      lower2(1,0) = 3;
      lower2(2,0) = 4;
      lower2(2,1) = 5;

      swap( lower1, lower2 );

      checkRows    ( lower1, 3UL );
      checkColumns ( lower1, 3UL );
      checkCapacity( lower1, 6UL );
      checkNonZeros( lower1, 6UL );
      checkNonZeros( lower1, 0UL, 1UL );
      checkNonZeros( lower1, 1UL, 2UL );
      checkNonZeros( lower1, 2UL, 3UL );

      if( lower1(0,0) != 1 || lower1(0,1) != 0 || lower1(0,2) != 0 ||
          lower1(1,0) != 3 || lower1(1,1) != 1 || lower1(1,2) != 0 ||
          lower1(2,0) != 4 || lower1(2,1) != 5 || lower1(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower1 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 3 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( lower2, 2UL );
      checkColumns ( lower2, 2UL );
      checkCapacity( lower2, 3UL );
      checkNonZeros( lower2, 3UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );

      if( lower2(0,0) != 1 || lower2(0,1) != 0 || lower2(1,0) != 2 || lower2(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix swap";

      OLT lower1( 2UL );
      lower1(1,0) = 2;

      OLT lower2( 3UL );
      lower2(1,0) = 3;
      lower2(2,0) = 4;
      lower2(2,1) = 5;

      swap( lower1, lower2 );

      checkRows    ( lower1, 3UL );
      checkColumns ( lower1, 3UL );
      checkCapacity( lower1, 6UL );
      checkNonZeros( lower1, 6UL );
      checkNonZeros( lower1, 0UL, 3UL );
      checkNonZeros( lower1, 1UL, 2UL );
      checkNonZeros( lower1, 2UL, 1UL );

      if( lower1(0,0) != 1 || lower1(0,1) != 0 || lower1(0,2) != 0 ||
          lower1(1,0) != 3 || lower1(1,1) != 1 || lower1(1,2) != 0 ||
          lower1(2,0) != 4 || lower1(2,1) != 5 || lower1(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower1 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 3 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( lower2, 2UL );
      checkColumns ( lower2, 2UL );
      checkCapacity( lower2, 3UL );
      checkNonZeros( lower2, 3UL );
      checkNonZeros( lower2, 0UL, 2UL );
      checkNonZeros( lower2, 1UL, 1UL );

      if( lower2(0,0) != 1 || lower2(0,1) != 0 || lower2(1,0) != 2 || lower2(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::find()";

      typedef LT::ConstIterator  ConstIterator;

      // Initialization check
      LT lower( 8UL, 10UL );
      lower(2,1) = 2;
      lower(4,2) = 3;

      checkRows    ( lower,  8UL );
      checkColumns ( lower,  8UL );
      checkCapacity( lower, 10UL );
      checkNonZeros( lower, 10UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 1UL );
      checkNonZeros( lower,  2UL, 2UL );
      checkNonZeros( lower,  3UL, 1UL );
      checkNonZeros( lower,  4UL, 2UL );
      checkNonZeros( lower,  5UL, 1UL );
      checkNonZeros( lower,  6UL, 1UL );
      checkNonZeros( lower,  7UL, 1UL );

      // Searching for the first element
      {
         ConstIterator pos( lower.find( 2UL, 1UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 1 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( lower.find( 4UL, 2UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (4,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a diagonal element
      {
         ConstIterator pos( lower.find( 6UL, 6UL ) );

         if( pos == lower.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,6)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 6 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 6\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( lower.find( 4UL, 0UL ) );

         if( pos != lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::find()";

      typedef OLT::ConstIterator  ConstIterator;

      // Initialization check
      OLT lower( 8UL, 10UL );
      lower(2,1) = 2;
      lower(4,2) = 3;

      checkRows    ( lower,  8UL );
      checkColumns ( lower,  8UL );
      checkCapacity( lower, 10UL );
      checkNonZeros( lower, 10UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 2UL );
      checkNonZeros( lower,  2UL, 2UL );
      checkNonZeros( lower,  3UL, 1UL );
      checkNonZeros( lower,  4UL, 1UL );
      checkNonZeros( lower,  5UL, 1UL );
      checkNonZeros( lower,  6UL, 1UL );
      checkNonZeros( lower,  7UL, 1UL );

      // Searching for the first element
      {
         ConstIterator pos( lower.find( 2UL, 1UL ) );

         if( pos == lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << lower << "\n";
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
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( lower.find( 4UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (4,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a diagonal element
      {
         ConstIterator pos( lower.find( 6UL, 6UL ) );

         if( pos == lower.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,6)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 6 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 6\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( lower.find( 4UL, 0UL ) );

         if( pos != lower.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::lowerBound()";

      typedef LT::ConstIterator  ConstIterator;

      // Initialization check
      LT lower( 6UL, 7UL );
      lower(4,2) = 2;

      checkRows    ( lower, 6UL );
      checkColumns ( lower, 6UL );
      checkCapacity( lower, 7UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );
      checkNonZeros( lower, 4UL, 2UL );
      checkNonZeros( lower, 5UL, 1UL );

      // Determining the lower bound for position (4,1)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 1UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << lower << "\n";
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
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,2)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 2UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,2)\n"
                << "   Current matrix:\n" << lower << "\n";
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
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,3)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 3UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,3)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,4)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 4UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,4)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,5)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 5UL ) );

         if( pos != lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,5)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::lowerBound()";

      typedef OLT::ConstIterator  ConstIterator;

      // Initialization check
      OLT lower( 6UL, 7UL );
      lower(4,2) = 2;

      checkRows    ( lower, 6UL );
      checkColumns ( lower, 6UL );
      checkCapacity( lower, 7UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );
      checkNonZeros( lower, 4UL, 1UL );
      checkNonZeros( lower, 5UL, 1UL );

      // Determining the lower bound for position (1,2)
      {
         ConstIterator pos( lower.lowerBound( 1UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,2)
      {
         ConstIterator pos( lower.lowerBound( 2UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (3,2)
      {
         ConstIterator pos( lower.lowerBound( 3UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,2)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (5,2)
      {
         ConstIterator pos( lower.lowerBound( 5UL, 2UL ) );

         if( pos != lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::upperBound()";

      typedef LT::ConstIterator  ConstIterator;

      // Initialization check
      LT lower( 6UL, 7UL );
      lower(4,2) = 2;

      checkRows    ( lower, 6UL );
      checkColumns ( lower, 6UL );
      checkCapacity( lower, 7UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );
      checkNonZeros( lower, 4UL, 2UL );
      checkNonZeros( lower, 5UL, 1UL );

      // Determining the upper bound for position (4,1)
      {
         ConstIterator pos( lower.upperBound( 4UL, 1UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << lower << "\n";
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
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,2)
      {
         ConstIterator pos( lower.upperBound( 4UL, 2UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,3)
      {
         ConstIterator pos( lower.upperBound( 4UL, 3UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,3)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,4)
      {
         ConstIterator pos( lower.upperBound( 4UL, 4UL ) );

         if( pos != lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,4)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,5)
      {
         ConstIterator pos( lower.upperBound( 4UL, 5UL ) );

         if( pos != lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,5)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::upperBound()";

      typedef OLT::ConstIterator  ConstIterator;

      // Initialization check
      OLT lower( 6UL, 7UL );
      lower(4,2) = 2;

      checkRows    ( lower, 6UL );
      checkColumns ( lower, 6UL );
      checkCapacity( lower, 7UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );
      checkNonZeros( lower, 4UL, 1UL );
      checkNonZeros( lower, 5UL, 1UL );

      // Determining the upper bound for position (1,2)
      {
         ConstIterator pos( lower.upperBound( 1UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,2)
      {
         ConstIterator pos( lower.upperBound( 2UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (3,2)
      {
         ConstIterator pos( lower.upperBound( 3UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,2)
      {
         ConstIterator pos( lower.upperBound( 4UL, 2UL ) );

         if( pos != lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (5,2)
      {
         ConstIterator pos( lower.upperBound( 5UL, 2UL ) );

         if( pos != lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         LT lower;

         if( isDefault( lower ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         LT lower( 3UL );

         if( isDefault( lower(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         LT lower( 3UL );
         lower(1,0) = 1;

         if( isDefault( lower(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
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
         OLT lower;

         if( isDefault( lower ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         OLT lower( 3UL );

         if( isDefault( lower(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         OLT lower( 3UL );
         lower(1,0) = 1;

         if( isDefault( lower(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c submatrix() function with the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c submatrix() function with the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      typedef blaze::Submatrix<LT>  SMT;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      SMT sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      if( sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(1,1) << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,0) = -5;

      if( sm(0,0) !=  1 || sm(0,1) != 0 ||
          sm(1,0) != -5 || sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  1  0 )\n( -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4  1  0 )\n(  7 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != 1 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      typedef blaze::Submatrix<OLT>  SMT;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      SMT sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      if( sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(1,1) << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,0) = -5;

      if( sm(0,0) !=  1 || sm(0,1) != 0 ||
          sm(1,0) != -5 || sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  1  0 )\n( -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4  1  0 )\n(  7 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != 1 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c row() function with the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c row() function with the UniLowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testRow()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      typedef blaze::SparseRow<LT>  RT;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      RT row1 = row( lower, 1UL );

      if( row1[0] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[0] << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it = row1.begin();

      if( it == row1.end() || it->value() != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      row1[0] = -5;

      if( row1[0] != -5 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -5 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -5 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -5 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      typedef blaze::SparseRow<OLT>  RT;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      RT row1 = row( lower, 1UL );

      if( row1[0] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[0] << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it = row1.begin();

      if( it == row1.end() || it->value() != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      row1[0] = -5;

      if( row1[0] != -5 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -5 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -5 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -5 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c column() function with the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c column() function with the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testColumn()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      typedef blaze::SparseColumn<LT>  CT;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      CT col0 = column( lower, 0UL );

      if( col0[0] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col0[0] << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it = col0.begin();

      if( it == col0.end() || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      col0[1] = -5;

      if( col0[0] != 1 || col0[1] != -5 || col0[2] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -5 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -5 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col0 );

      if( col0[0] != 1 || col0[1] != 0 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      typedef blaze::SparseColumn<OLT>  CT;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      CT col0 = column( lower, 0UL );

      if( col0[0] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col0[0] << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it = col0.begin();

      if( it == col0.end() || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      col0[1] = -5;

      if( col0[0] != 1 || col0[1] != -5 || col0[2] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 1 -5  7 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -5 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -5 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col0 );

      if( col0[0] != 1 || col0[1] != 0 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace unilowermatrix

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
   std::cout << "   Running UniLowerMatrix sparse test..." << std::endl;

   try
   {
      RUN_UNILOWERMATRIX_SPARSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UniLowerMatrix sparse test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
