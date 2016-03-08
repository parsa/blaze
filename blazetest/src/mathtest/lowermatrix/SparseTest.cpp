//=================================================================================================
/*!
//  \file src/mathtest/lowermatrix/SparseTest.cpp
//  \brief Source file for the LowerMatrix sparse test
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
#include <blaze/math/SparseSubmatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/lowermatrix/SparseTest.h>


namespace blazetest {

namespace mathtest {

namespace lowermatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the LowerMatrix sparse test.
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
   testScaling();
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
/*!\brief Test of the LowerMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the LowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (CompressedMatrix)
   {
      test_ = "Row-major LowerMatrix default constructor (CompressedMatrix)";

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
      test_ = "Row-major LowerMatrix size constructor (CompressedMatrix)";

      const LT lower( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkNonZeros( lower, 0UL );
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Row-major LowerMatrix copy constructor (0x0)";

      const LT lower1;
      const LT lower2( lower1 );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Row-major LowerMatrix copy constructor (3x3)";

      LT lower1( 3UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      const LT lower2( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Row-major LowerMatrix move constructor (0x0)";

      LT lower1;
      LT lower2( std::move( lower1 ) );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Row-major LowerMatrix move constructor (3x3)";

      LT lower1( 3UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      LT lower2( std::move( lower1 ) );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Row-major LowerMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      const LT lower( mat );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Conversion constructor (lower)
   {
      test_ = "Row-major LowerMatrix conversion constructor (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      const LT lower( mat );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-lower)
   {
      test_ = "Row-major LowerMatrix conversion constructor (non-lower)";

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
             << " Error: Setup of non-lower LowerMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (LowerMatrix)
   {
      test_ = "Row-major LowerMatrix conversion constructor (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      const LT lower2( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (CompressedMatrix)
   {
      test_ = "Column-major LowerMatrix default constructor (CompressedMatrix)";

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
      test_ = "Column-major LowerMatrix size constructor (CompressedMatrix)";

      const OLT lower( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkNonZeros( lower, 0UL );
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Column-major LowerMatrix copy constructor (0x0)";

      const OLT lower1;
      const OLT lower2( lower1 );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Column-major LowerMatrix copy constructor (3x3)";

      OLT lower1( 3UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      const OLT lower2( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Column-major LowerMatrix move constructor (0x0)";

      OLT lower1;
      OLT lower2( std::move( lower1 ) );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Column-major LowerMatrix move constructor (3x3)";

      OLT lower1( 3UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      OLT lower2( std::move( lower1 ) );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Column-major LowerMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      const OLT lower( mat );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Conversion constructor (lower)
   {
      test_ = "Column-major LowerMatrix conversion constructor (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      const OLT lower( mat );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-lower)
   {
      test_ = "Column-major LowerMatrix conversion constructor (non-lower)";

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
             << " Error: Setup of non-lower LowerMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (LowerMatrix)
   {
      test_ = "Column-major LowerMatrix conversion constructor (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      const OLT lower2( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LowerMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the LowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAssignment()
{
   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Row-major LowerMatrix copy assignment (0x0)";

      LT lower1, lower2;

      lower2 = lower1;

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Row-major LowerMatrix copy assignment (3x3)";

      LT lower1( 3UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Row-major LowerMatrix move assignment (0x0)";

      LT lower1, lower2;

      lower2 = std::move( lower1 );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Row-major LowerMatrix move assignment (3x3)";

      LT lower1( 3UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      LT lower2;
      lower2 = std::move( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major LowerMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      LT lower;
      lower = mat;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Row-major/row-major dense matrix assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      LT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      LT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix assignment (non-lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         LT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (non-lower)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix assignment (non-lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         LT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major LowerMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      LT lower;
      lower = mat;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
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
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
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
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         LT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (non-lower)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         LT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::CompressedMatrix<unsigned int,blaze::rowMajor> > lower1( 3UL, 5UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > lower1( 3UL, 5UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      LT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Column-major LowerMatrix copy assignment (0x0)";

      OLT lower1, lower2;

      lower2 = lower1;

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Column-major LowerMatrix copy assignment (3x3)";

      OLT lower1( 3UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Column-major LowerMatrix move assignment (0x0)";

      OLT lower1, lower2;

      lower2 = std::move( lower1 );

      checkRows    ( lower2, 0UL );
      checkColumns ( lower2, 0UL );
      checkNonZeros( lower2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Column-major LowerMatrix move assignment (3x3)";

      OLT lower1( 3UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      OLT lower2;
      lower2 = std::move( lower1 );

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major LowerMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      OLT lower;
      lower = mat;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Column-major/row-major dense matrix assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      OLT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      OLT lower;
      lower = mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix assignment (non-lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         OLT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (non-lower)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix assignment (non-lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         OLT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major LowerMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat;

      OLT lower;
      lower = mat;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
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
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
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
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         OLT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (non-lower)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) =  1;
      mat(0,2) =  5;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;

      try {
         OLT lower;
         lower = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > lower1( 3UL, 5UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::CompressedMatrix<unsigned int,blaze::columnMajor> > lower1( 3UL, 5UL );
      lower1(0,0) =  1;
      lower1(1,0) = -4;
      lower1(1,1) =  2;
      lower1(2,0) =  7;
      lower1(2,2) =  3;

      OLT lower2;
      lower2 = lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -4 || lower2(1,1) != 2 || lower2(1,2) != 0 ||
          lower2(2,0) !=  7 || lower2(2,1) != 0 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 2 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LowerMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   // Row-major/row-major dense matrix addition assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) =  2;
      mat(1,1) = -2;
      mat(2,0) =  6;
      mat(2,1) =  5;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 13 || lower(2,1) != 5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) =  2;
      mat(1,1) = -2;
      mat(2,0) =  6;
      mat(2,1) =  5;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 13 || lower(2,1) != 5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix addition assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix addition assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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
      test_ = "Row-major/column-major LowerMatrix dense matrix addition assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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

   // Row-major/row-major dense matrix addition assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix addition assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(1,0) =  2;
      lower1(1,1) = -2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 += lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) != 0 || lower2(1,2) != 0 ||
          lower2(2,0) != 13 || lower2(2,1) != 5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix addition assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(1,0) =  2;
      lower1(1,1) = -2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 += lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) != 0 || lower2(1,2) != 0 ||
          lower2(2,0) != 13 || lower2(2,1) != 5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix addition assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(1,0) =  2;
      mat(1,1) = -2;
      mat(2,0) =  6;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 13 || lower(2,1) != 5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(1,0) =  2;
      mat(1,1) = -2;
      mat(2,0) =  6;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 13 || lower(2,1) != 5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix addition assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix addition assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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
      test_ = "Row-major/column-major LowerMatrix sparse matrix addition assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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

   // Row-major/row-major sparse matrix addition assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix addition assignment (LowerMatrix)";

      LT lower1( 3UL, 4UL );
      lower1(1,0) =  2;
      lower1(1,1) = -2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 += lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) != 0 || lower2(1,2) != 0 ||
          lower2(2,0) != 13 || lower2(2,1) != 5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix addition assignment (LowerMatrix)";

      OLT lower1( 3UL, 4UL );
      lower1(1,0) =  2;
      lower1(1,1) = -2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 += lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) != 0 || lower2(1,2) != 0 ||
          lower2(2,0) != 13 || lower2(2,1) != 5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   // Column-major/row-major dense matrix addition assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) =  2;
      mat(1,1) = -2;
      mat(2,0) =  6;
      mat(2,1) =  5;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 13 || lower(2,1) != 5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) =  2;
      mat(1,1) = -2;
      mat(2,0) =  6;
      mat(2,1) =  5;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 13 || lower(2,1) != 5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix addition assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix addition assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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
      test_ = "Column-major/column-major LowerMatrix dense matrix addition assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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

   // Column-major/row-major dense matrix addition assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix addition assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(1,0) =  2;
      lower1(1,1) = -2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 += lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) != 0 || lower2(1,2) != 0 ||
          lower2(2,0) != 13 || lower2(2,1) != 5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix addition assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(1,0) =  2;
      lower1(1,1) = -2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 += lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) != 0 || lower2(1,2) != 0 ||
          lower2(2,0) != 13 || lower2(2,1) != 5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix addition assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(1,0) =  2;
      mat(1,1) = -2;
      mat(2,0) =  6;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 13 || lower(2,1) != 5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(1,0) =  2;
      mat(1,1) = -2;
      mat(2,0) =  6;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower += mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 13 || lower(2,1) != 5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix addition assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix addition assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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
      test_ = "Column-major/column-major LowerMatrix sparse matrix addition assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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

   // Column-major/row-major sparse matrix addition assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix addition assignment (LowerMatrix)";

      LT lower1( 3UL, 4UL );
      lower1(1,0) =  2;
      lower1(1,1) = -2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 += lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) != 0 || lower2(1,2) != 0 ||
          lower2(2,0) != 13 || lower2(2,1) != 5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix addition assignment (LowerMatrix)";

      OLT lower1( 3UL, 4UL );
      lower1(1,0) =  2;
      lower1(1,1) = -2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 += lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) != 0 || lower2(1,2) != 0 ||
          lower2(2,0) != 13 || lower2(2,1) != 5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 0 0 )\n( 13 5 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LowerMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major dense matrix subtraction assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) = -2;
      mat(1,1) =  2;
      mat(2,0) =  6;
      mat(2,1) =  5;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  0 || lower(1,2) != 0 ||
          lower(2,0) !=  1 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) = -2;
      mat(1,1) =  2;
      mat(2,0) =  6;
      mat(2,1) =  5;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  0 || lower(1,2) != 0 ||
          lower(2,0) !=  1 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix subtraction assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix subtraction assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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
      test_ = "Row-major/column-major LowerMatrix dense matrix subtraction assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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

   // Row-major/row-major dense matrix subtraction assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix subtraction assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(1,0) = -2;
      lower1(1,1) =  2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 -= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) !=  0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) !=  0 || lower2(1,2) != 0 ||
          lower2(2,0) !=  1 || lower2(2,1) != -5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix subtraction assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(1,0) = -2;
      lower1(1,1) =  2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 -= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) !=  0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) !=  0 || lower2(1,2) != 0 ||
          lower2(2,0) !=  1 || lower2(2,1) != -5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(1,0) = -2;
      mat(1,1) =  2;
      mat(2,0) =  6;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  0 || lower(1,2) != 0 ||
          lower(2,0) !=  1 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(1,0) = -2;
      mat(1,1) =  2;
      mat(2,0) =  6;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  0 || lower(1,2) != 0 ||
          lower(2,0) !=  1 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix subtraction assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix subtraction assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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
      test_ = "Row-major/column-major LowerMatrix sparse matrix subtraction assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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

   // Row-major/row-major sparse matrix subtraction assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix subtraction assignment (LowerMatrix)";

      LT lower1( 3UL, 4UL );
      lower1(1,0) = -2;
      lower1(1,1) =  2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 -= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) !=  0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) !=  0 || lower2(1,2) != 0 ||
          lower2(2,0) !=  1 || lower2(2,1) != -5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix subtraction assignment (LowerMatrix)";

      OLT lower1( 3UL, 4UL );
      lower1(1,0) = -2;
      lower1(1,1) =  2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 -= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 3UL );

      if( lower2(0,0) !=  1 || lower2(0,1) !=  0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) !=  0 || lower2(1,2) != 0 ||
          lower2(2,0) !=  1 || lower2(2,1) != -5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major dense matrix subtraction assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) = -2;
      mat(1,1) =  2;
      mat(2,0) =  6;
      mat(2,1) =  5;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  0 || lower(1,2) != 0 ||
          lower(2,0) !=  1 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,0) = -2;
      mat(1,1) =  2;
      mat(2,0) =  6;
      mat(2,1) =  5;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  0 || lower(1,2) != 0 ||
          lower(2,0) !=  1 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix subtraction assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix subtraction assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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
      test_ = "Column-major/column-major LowerMatrix dense matrix subtraction assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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

   // Column-major/row-major dense matrix subtraction assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix subtraction assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(1,0) = -2;
      lower1(1,1) =  2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 -= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) !=  0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) !=  0 || lower2(1,2) != 0 ||
          lower2(2,0) !=  1 || lower2(2,1) != -5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix subtraction assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(1,0) = -2;
      lower1(1,1) =  2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 -= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) !=  0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) !=  0 || lower2(1,2) != 0 ||
          lower2(2,0) !=  1 || lower2(2,1) != -5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(1,0) = -2;
      mat(1,1) =  2;
      mat(2,0) =  6;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  0 || lower(1,2) != 0 ||
          lower(2,0) !=  1 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(1,0) = -2;
      mat(1,1) =  2;
      mat(2,0) =  6;
      mat(2,1) =  5;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower -= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) !=  0 || lower(1,2) != 0 ||
          lower(2,0) !=  1 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix subtraction assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix subtraction assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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
      test_ = "Column-major/column-major LowerMatrix sparse matrix subtraction assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

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

   // Column-major/row-major sparse matrix subtraction assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix subtraction assignment (LowerMatrix)";

      LT lower1( 3UL, 4UL );
      lower1(1,0) = -2;
      lower1(1,1) =  2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 -= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) !=  0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) !=  0 || lower2(1,2) != 0 ||
          lower2(2,0) !=  1 || lower2(2,1) != -5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix subtraction assignment (LowerMatrix)";

      OLT lower1( 3UL, 4UL );
      lower1(1,0) = -2;
      lower1(1,1) =  2;
      lower1(2,0) =  6;
      lower1(2,1) =  5;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 -= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 6UL );
      checkNonZeros( lower2, 6UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  1 || lower2(0,1) !=  0 || lower2(0,2) != 0 ||
          lower2(1,0) != -2 || lower2(1,1) !=  0 || lower2(1,2) != 0 ||
          lower2(2,0) !=  1 || lower2(2,1) != -5 || lower2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -2  0  0 )\n(  1 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LowerMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major dense matrix multiplication assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix multiplication assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix multiplication assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (non-lower)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix multiplication assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix multiplication assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix multiplication assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix multiplication assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix multiplication assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (non-lower)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix multiplication assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix multiplication assignment (LowerMatrix)";

      LT lower1( 3UL, 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix multiplication assignment (LowerMatrix)";

      OLT lower1( 3UL, 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major dense matrix multiplication assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix multiplication assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix multiplication assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (non-lower)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix multiplication assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix multiplication assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower1;
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix multiplication assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower1;
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 5UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix multiplication assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix multiplication assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (non-lower)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix multiplication assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix multiplication assignment (LowerMatrix)";

      LT lower1( 3UL, 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix multiplication assignment (LowerMatrix)";

      OLT lower1( 3UL, 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 5UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all LowerMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      LT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower *= 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      LT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower = lower * 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      LT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower = 2 * lower;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      LT lower( 3UL );
      lower(1,0) =  2;
      lower(2,0) = -4;
      lower(2,1) =  6;
      lower(2,2) = -8;

      lower /= 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      LT lower( 3UL );
      lower(1,0) =  2;
      lower(2,0) = -4;
      lower(2,1) =  6;
      lower(2,2) = -8;

      lower = lower / 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major LowerMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::scale()";

      // Initialization check
      LT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      lower.scale( 2 );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      lower.scale( 0.5 );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major LowerMatrix::scale() (complex)";

      using blaze::complex;

      blaze::LowerMatrix< blaze::CompressedMatrix<complex<float>,blaze::rowMajor> > lower( 2UL );
      lower(0,0) = complex<float>( 1.0F, 0.0F );
      lower(1,0) = complex<float>( 2.0F, 0.0F );
      lower(1,1) = complex<float>( 4.0F, 0.0F );

      lower.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );

      if( lower(0,0) != complex<float>( 3.0F, 0.0F ) || lower(0,1) != complex<float>(  0.0F, 0.0F ) ||
          lower(1,0) != complex<float>( 6.0F, 0.0F ) || lower(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( ( 3,0) ( 0,0)\n( 6,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      OLT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower *= 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      OLT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower = lower * 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      OLT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower = 2 * lower;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      OLT lower( 3UL );
      lower(1,0) =  2;
      lower(2,0) = -4;
      lower(2,1) =  6;
      lower(2,2) = -8;

      lower /= 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      OLT lower( 3UL );
      lower(1,0) =  2;
      lower(2,0) = -4;
      lower(2,1) =  6;
      lower(2,2) = -8;

      lower = lower / 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major LowerMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::scale()";

      // Initialization check
      OLT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      lower.scale( 2 );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      lower.scale( 0.5 );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major LowerMatrix::scale() (complex)";

      using blaze::complex;

      blaze::LowerMatrix< blaze::CompressedMatrix<complex<float>,blaze::columnMajor> > lower( 2UL );
      lower(0,0) = complex<float>( 1.0F, 0.0F );
      lower(1,0) = complex<float>( 2.0F, 0.0F );
      lower(1,1) = complex<float>( 4.0F, 0.0F );

      lower.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );

      if( lower(0,0) != complex<float>( 3.0F, 0.0F ) || lower(0,1) != complex<float>(  0.0F, 0.0F ) ||
          lower(1,0) != complex<float>( 6.0F, 0.0F ) || lower(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( ( 3,0) ( 0,0)\n( 6,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LowerMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the LowerMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::operator()";

      // Good cases
      {
         LT lower( 3UL );

         // Writing the diagonal element (1,1)
         lower(1,1) = 1;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (2,1)
         lower(2,1) = 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 2UL );
         checkNonZeros( lower, 2UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (1,0)
         lower(1,0) = lower(2,1);

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the lower element (2,0)
         lower(2,0) += 3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the lower element (1,0)
         lower(1,0) -= 4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0 0 0 )\n( -2 1 0 )\n(  3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the lower element (2,1)
         lower(2,1) *= -3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) !=  0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -6 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0  0  0 )\n( -2  1  0 )\n(  3 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the lower element (2,1)
         lower(2,1) /= 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) !=  0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -3 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0  0  0 )\n( -2  1  0 )\n(  3 -3  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         LT lower( 3UL );

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
      test_ = "Column-major LowerMatrix::operator()";

      // Good cases
      {
         OLT lower( 3UL );

         // Writing the diagonal element (1,1)
         lower(1,1) = 1;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (2,1)
         lower(2,1) = 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 2UL );
         checkNonZeros( lower, 2UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (1,0)
         lower(1,0) = lower(2,1);

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the lower element (2,0)
         lower(2,0) += 3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the lower element (1,0)
         lower(1,0) -= 4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0 0 0 )\n( -2 1 0 )\n(  3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the lower element (2,1)
         lower(2,1) *= -3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) !=  0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -6 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0  0  0 )\n( -2  1  0 )\n(  3 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the lower element (2,1)
         lower(2,1) /= 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) !=  0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -3 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0  0  0 )\n( -2  1  0 )\n(  3 -3  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         OLT lower( 3UL );

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
/*!\brief Test of the LowerMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the LowerMatrix
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
      lower(0,0) =  1;
      lower(1,1) = -2;
      lower(2,0) =  3;
      lower(2,2) =  4;

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

         if( it == end( lower, 1UL ) || it->value() != -2 ) {
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

         if( it == end || it->value() != 4 ) {
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

         int value = 7;

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it = value++;
         }

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 7 || lower(2,1) !=  0 || lower(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 7  0  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Row-major addition assignment to lower elements via Iterator";

         int value = 4;

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it += value++;
         }

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 ||
             lower(1,0) !=  0 || lower(1,1) != -2 || lower(1,2) !=  0 ||
             lower(2,0) != 11 || lower(2,1) !=  0 || lower(2,2) != 13 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n(  0 -2  0 )\n( 11  0 13 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Row-major subtraction assignment to lower elements via Iterator";

         int value = 4;

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it -= value++;
         }

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 7 || lower(2,1) !=  0 || lower(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 7  0  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Row-major multiplication assignment to lower elements via Iterator";

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it *= 2;
         }

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 ||
             lower(1,0) !=  0 || lower(1,1) != -2 || lower(1,2) !=  0 ||
             lower(2,0) != 14 || lower(2,1) !=  0 || lower(2,2) != 16 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n(  0 -2  0 )\n( 14  0 16 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Row-major division assignment to lower elements via Iterator";

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it /= 2;
         }

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 7 || lower(2,1) !=  0 || lower(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 7  0  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      typedef OLT::Iterator       Iterator;
      typedef OLT::ConstIterator  ConstIterator;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,1) = -2;
      lower(2,0) =  3;
      lower(2,2) =  4;

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

         if( it == end( lower, 1UL ) || it->value() != -2 ) {
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

         int value = 7;

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it = value++;
         }

         if( lower(0,0) != 7 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 8 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 7  0  0 )\n( 0 -2  0 )\n( 8  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Column-major addition assignment to lower elements via Iterator";

         int value = 4;

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it += value++;
         }

         if( lower(0,0) != 11 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 13 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 11  0  0 )\n(  0 -2  0 )\n( 13  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Column-major subtraction assignment to lower elements via Iterator";

         int value = 4;

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it -= value++;
         }

         if( lower(0,0) != 7 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 8 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 7  0  0 )\n( 0 -2  0 )\n( 8  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Column-major multiplication assignment to lower elements via Iterator";

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it *= 2;
         }

         if( lower(0,0) != 14 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 16 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 14  0  0 )\n(  0 -2  0 )\n( 16  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Column-major division assignment to lower elements via Iterator";

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it /= 2;
         }

         if( lower(0,0) != 7 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 8 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 7  0  0 )\n( 0 -2  0 )\n( 8  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::nonZeros()";

      // Empty matrix
      {
         LT lower( 3UL );

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkNonZeros( lower, 0UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         LT lower( 3UL );
         lower(0,0) =  1;
         lower(1,1) = -2;
         lower(2,1) =  3;
         lower(2,2) = -4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) !=  0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) !=  0 ||
             lower(2,0) != 0 || lower(2,1) !=  3 || lower(2,2) != -4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  3 -4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         LT lower( 3UL );
         lower(0,0) = -1;
         lower(1,0) =  2;
         lower(1,1) =  3;
         lower(2,0) = -4;
         lower(2,1) = -5;
         lower(2,2) =  6;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) != -1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  2 || lower(1,1) !=  3 || lower(1,2) != 0 ||
             lower(2,0) != -4 || lower(2,1) != -5 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( -1  0  0 )\n(  2  3  0 )\n( -4 -5  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::nonZeros()";

      // Empty matrix
      {
         OLT lower( 3UL );

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkNonZeros( lower, 0UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         OLT lower( 3UL );
         lower(0,0) =  1;
         lower(1,1) = -2;
         lower(2,1) =  3;
         lower(2,2) = -4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 4UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) !=  0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) !=  0 ||
             lower(2,0) != 0 || lower(2,1) !=  3 || lower(2,2) != -4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  3 -4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         OLT lower( 3UL );
         lower(0,0) = -1;
         lower(1,0) =  2;
         lower(1,1) =  3;
         lower(2,0) = -4;
         lower(2,1) = -5;
         lower(2,2) =  6;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 6UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != -1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  2 || lower(1,1) !=  3 || lower(1,2) != 0 ||
             lower(2,0) != -4 || lower(2,1) != -5 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( -1  0  0 )\n(  2  3  0 )\n( -4 -5  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::reset()";

      // Initialization check
      LT lower( 3UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(2,2) = 6;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( lower(1,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( lower(0,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( lower, 1UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );

      if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::reset()";

      // Initialization check
      OLT lower( 3UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(2,2) = 6;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( lower(1,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( lower(0,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting column 1
      reset( lower, 1UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );

      if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::clear()";

      // Initialization check
      LT lower( 3UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(2,2) = 6;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(1,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
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
      test_ = "Column-major LowerMatrix::clear()";

      // Initialization check
      OLT lower( 3UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(2,2) = 6;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(1,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 6UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
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
/*!\brief Test of the \c set() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::set()";

      typedef LT::Iterator  Iterator;

      // Initialization check
      LT lower( 4UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = lower.set( 2UL, 1UL, 1 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( pos->value() != 1 || pos->index() != 1UL ) {
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

         if( lower(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = lower.set( 2UL, 2UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 2UL );
         checkNonZeros( lower, 2UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 0UL );

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

         if( lower(2,1) != 1 || lower(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = lower.set( 2UL, 0UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 3UL );
         checkNonZeros( lower, 3UL, 0UL );

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

         if( lower(2,0) != 3 || lower(2,1) != 1 || lower(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = lower.set( 2UL, 1UL, 4 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 3UL );
         checkNonZeros( lower, 3UL, 0UL );

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

         if( lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 3 4 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::set()";

      typedef OLT::Iterator  Iterator;

      // Initialization check
      OLT lower( 4UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = lower.set( 2UL, 1UL, 1 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( pos->value() != 1 || pos->index() != 2UL ) {
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

         if( lower(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = lower.set( 1UL, 1UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 2UL );
         checkNonZeros( lower, 2UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(1,1) != 2 || lower(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = lower.set( 3UL, 1UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

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

         if( lower(1,1) != 2 || lower(2,1) != 1 || lower(3,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 1 0 0 )\n( 0 3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = lower.set( 2UL, 1UL, 4 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

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

         if( lower(1,1) != 2 || lower(2,1) != 4 || lower(3,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 4 0 0 )\n( 0 3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::insert()";

      typedef LT::Iterator  Iterator;

      // Initialization check
      LT lower( 4UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = lower.insert( 2UL, 1UL, 1 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( pos->value() != 1 || pos->index() != 1UL ) {
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

         if( lower(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = lower.insert( 2UL, 2UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 2UL );
         checkNonZeros( lower, 2UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 0UL );

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

         if( lower(2,1) != 1 || lower(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = lower.insert( 2UL, 0UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 3UL );
         checkNonZeros( lower, 3UL, 0UL );

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

         if( lower(2,0) != 3 || lower(2,1) != 1 || lower(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
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
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::insert()";

      typedef OLT::Iterator  Iterator;

      // Initialization check
      OLT lower( 4UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = lower.insert( 2UL, 1UL, 1 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( pos->value() != 1 || pos->index() != 2UL ) {
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

         if( lower(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = lower.insert( 1UL, 1UL, 2 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 2UL );
         checkNonZeros( lower, 2UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(1,1) != 2 || lower(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = lower.insert( 3UL, 1UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 3UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

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

         if( lower(1,1) != 2 || lower(2,1) != 1 || lower(3,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 1 0 0 )\n( 0 3 0 0 )\n";
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
             << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 1 0 0 )\n( 0 3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::append()";

      // Appending with pre-allocation in each row
      {
         // Initialization check
         LT lower( 4UL, 5UL );
         lower.reserve( 0UL, 1UL );
         lower.reserve( 2UL, 2UL );
         lower.reserve( 3UL, 2UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkNonZeros( lower, 0UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         // Appending one non-zero element
         lower.append( 2UL, 1UL, 1 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( lower(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 0UL, 0UL, 2 );
         lower.append( 3UL, 0UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 2 || lower(2,1) != 1 || lower(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 3UL, 2UL, 4 );
         lower.append( 2UL, 2UL, 5 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 2UL );

         if( lower(0,0) != 2 || lower(2,1) != 1 || lower(2,2) != 5 || lower(3,0) != 3 || lower(3,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 0 0 )\n( 0 1 5 0 )\n( 3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         // Initialization check
         LT lower( 4UL, 5UL );
         lower.reserve( 0UL, 1UL );
         lower.reserve( 2UL, 2UL );
         lower.reserve( 3UL, 2UL );

         // Appending one non-zero element
         lower.append( 0UL, 0UL, 1 );
         lower.finalize( 0UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( lower(0,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 2UL, 1UL, 2 );
         lower.append( 2UL, 2UL, 3 );
         lower.finalize( 2UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( lower(0,0) != 1 || lower(2,1) != 2 || lower(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 2 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 3UL, 0UL, 4 );
         lower.append( 3UL, 2UL, 5 );
         lower.finalize( 3UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 2UL );

         if( lower(0,0) != 1 || lower(2,1) != 2 || lower(2,2) != 3 || lower(3,0) != 4 || lower(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 2 3 0 )\n( 4 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::append()";

      // Appending with pre-allocation in each column
      {
         // Initialization check
         OLT lower( 4UL, 5UL );
         lower.reserve( 0UL, 2UL );
         lower.reserve( 1UL, 2UL );
         lower.reserve( 3UL, 1UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkNonZeros( lower, 0UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         // Appending one non-zero element
         lower.append( 1UL, 1UL, 1 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( lower(1,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 1UL, 0UL, 2 );
         lower.append( 3UL, 3UL, 3 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(1,0) != 2 || lower(1,1) != 1 || lower(3,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 2 1 0 0 )\n( 0 0 0 0 )\n( 0 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 3UL, 0UL, 4 );
         lower.append( 2UL, 1UL, 5 );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(1,0) != 2 || lower(1,1) != 1 || lower(2,1) != 5 || lower(3,0) != 4 || lower(3,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 2 1 0 0 )\n( 0 5 0 0 )\n( 4 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with column finalization
      {
         // Initialization check
         OLT lower( 4UL, 5UL );
         lower.reserve( 0UL, 1UL );
         lower.reserve( 1UL, 2UL );
         lower.reserve( 2UL, 2UL );

         // Appending one non-zero element
         lower.append( 1UL, 0UL, 1 );
         lower.finalize( 0UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 1UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( lower(1,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 1UL, 1UL, 2 );
         lower.append( 3UL, 1UL, 3 );
         lower.finalize( 1UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 3UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( lower(1,0) != 1 || lower(1,1) != 2 || lower(3,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 2 0 0 )\n( 0 0 0 0 )\n( 0 3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         lower.append( 2UL, 2UL, 4 );
         lower.append( 3UL, 2UL, 5 );
         lower.finalize( 2UL );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 5UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 0UL );

         if( lower(1,0) != 1 || lower(1,1) != 2 || lower(2,2) != 4 || lower(3,1) != 3 || lower(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 2 0 0 )\n( 0 0 4 0 )\n( 0 3 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::erase( size_t, size_t )";

      // Initialization check
      LT lower( 4UL, 8UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(3,0) = 6;
      lower(3,1) = 7;
      lower(3,3) = 8;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 3 ||
          lower(2,0) != 4 || lower(2,1) != 5 ||
          lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 3 0 0 )\n( 4 5 0 0 )\n( 6 7 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,0)
      lower.erase( 1UL, size_t(0) );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,1) != 3 ||
          lower(2,0) != 4 || lower(2,1) != 5 ||
          lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 5 0 0 )\n( 6 7 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      lower.erase( 2UL, 1UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,1) != 3 ||
          lower(2,0) != 4 ||
          lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 6 7 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (3,1)
      lower.erase( 3UL, 1UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 2UL );

      if( lower(0,0) != 1 ||
          lower(1,1) != 3 ||
          lower(2,0) != 4 ||
          lower(3,0) != 6 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 6 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      lower.erase( 3UL, 2UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 2UL );

      if( lower(0,0) != 1 ||
          lower(1,1) != 3 ||
          lower(2,0) != 4 ||
          lower(3,0) != 6 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 6 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::erase( size_t, Iterator )";

      typedef LT::Iterator  Iterator;

      // Initialization check
      LT lower( 4UL, 8UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(3,0) = 6;
      lower(3,1) = 7;
      lower(3,3) = 8;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 3 ||
          lower(2,0) != 4 || lower(2,1) != 5 ||
          lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 3 0 0 )\n( 4 5 0 0 )\n( 6 7 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,0)
      {
         Iterator pos = lower.erase( 1UL, lower.find( 1UL, 0UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 3UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 3 ||
             lower(2,0) != 4 || lower(2,1) != 5 ||
             lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 5 0 0 )\n( 6 7 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 3 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (2,1)
      {
         Iterator pos = lower.erase( 2UL, lower.find( 2UL, 1UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 3UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 3 ||
             lower(2,0) != 4 ||
             lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 6 7 0 8 )\n";
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

      // Erasing the element at (3,1)
      {
         Iterator pos = lower.erase( 3UL, lower.find( 3UL, 1UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 2UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 3 ||
             lower(2,0) != 4 ||
             lower(3,0) != 6 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 6 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 8 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 8\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = lower.erase( 3UL, lower.find( 3UL, 2UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 2UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 3 ||
             lower(2,0) != 4 ||
             lower(3,0) != 6 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 6 0 0 8 )\n";
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
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::erase( size_t, Iterator, Iterator )";

      typedef LT::Iterator  Iterator;

      // Initialization check
      LT lower( 4UL, 8UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(3,0) = 6;
      lower(3,1) = 7;
      lower(3,3) = 8;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( lower(0,0) != 1 ||
          lower(1,0) != 2 || lower(1,1) != 3 ||
          lower(2,0) != 4 || lower(2,1) != 5 ||
          lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 3 0 0 )\n( 4 5 0 0 )\n( 6 7 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from the beginning of row 1 to (1,1)
      {
         Iterator pos = lower.erase( 1UL, lower.begin( 1UL ), lower.find( 1UL, 1UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 3UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 3 ||
             lower(2,0) != 4 || lower(2,1) != 5 ||
             lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 5 0 0 )\n( 6 7 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 3 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (2,1) to the row end
      {
         Iterator pos = lower.erase( 2UL, lower.find( 2UL, 1UL ), lower.end( 2UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 3UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 3 ||
             lower(2,0) != 4 ||
             lower(3,0) != 6 || lower(3,1) != 7 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 6 7 0 8 )\n";
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

      // Erasing the elements from (3,0) to (3,3)
      {
         Iterator pos = lower.erase( 3UL, lower.find( 3UL, 0UL ), lower.find( 3UL, 3UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 3 ||
             lower(2,0) != 4 ||
             lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a multi-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 8 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 8\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = lower.erase( 3UL, lower.find( 3UL, 3UL ), lower.find( 3UL, 3UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(1,1) != 3 ||
             lower(2,0) != 4 ||
             lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 3 0 0 )\n( 4 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 8 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 8\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::erase( size_t, size_t )";

      // Initialization check
      OLT lower( 4UL, 8UL );
      lower(0,0) = 1;
      lower(2,0) = 2;
      lower(2,1) = 3;
      lower(2,2) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;
      lower(3,2) = 7;
      lower(3,3) = 8;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(2,0) != 2 || lower(2,1) != 3 || lower(2,2) != 4 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,2) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 3 4 0 )\n( 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      lower.erase( 2UL, 1UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 7UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(2,0) != 2 || lower(2,2) != 4 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,2) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 0 4 0 )\n( 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (3,2)
      lower.erase( 3UL, 2UL );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(2,0) != 2 || lower(2,2) != 4 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 0 4 0 )\n( 5 6 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,0)
      lower.erase( 2UL, size_t(0) );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(2,2) != 4 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 4 0 )\n( 5 6 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      lower.erase( 1UL, size_t(0) );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(2,2) != 4 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 4 0 )\n( 5 6 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::erase( size_t, Iterator )";

      typedef OLT::Iterator  Iterator;

      // Initialization check
      OLT lower( 4UL, 8UL );
      lower(0,0) = 1;
      lower(2,0) = 2;
      lower(2,1) = 3;
      lower(2,2) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;
      lower(3,2) = 7;
      lower(3,3) = 8;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(2,0) != 2 || lower(2,1) != 3 || lower(2,2) != 4 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,2) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 3 4 0 )\n( 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      {
         Iterator pos = lower.erase( 1UL, lower.find( 2UL, 1UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(2,0) != 2 || lower(2,2) != 4 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,2) != 7 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 0 4 0 )\n( 5 6 7 8 )\n";
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

      // Erasing the element at (3,2)
      {
         Iterator pos = lower.erase( 2UL, lower.find( 3UL, 2UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(2,0) != 2 || lower(2,2) != 4 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 0 4 0 )\n( 5 6 0 8 )\n";
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

      // Erasing the element at (2,0)
      {
         Iterator pos = lower.erase( 0UL, lower.find( 2UL, 0UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(2,2) != 4 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 4 0 )\n( 5 6 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = lower.erase( 0UL, lower.find( 1UL, 0UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(2,2) != 4 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 4 0 )\n( 5 6 0 8 )\n";
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
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::erase( size_t, Iterator, Iterator )";

      typedef OLT::Iterator  Iterator;

      // Initialization check
      OLT lower( 4UL, 8UL );
      lower(0,0) = 1;
      lower(2,0) = 2;
      lower(2,1) = 3;
      lower(2,2) = 4;
      lower(3,0) = 5;
      lower(3,1) = 6;
      lower(3,2) = 7;
      lower(3,3) = 8;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 8UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( lower(0,0) != 1 ||
          lower(2,0) != 2 || lower(2,1) != 3 || lower(2,2) != 4 ||
          lower(3,0) != 5 || lower(3,1) != 6 || lower(3,2) != 7 || lower(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 3 4 0 )\n( 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from the beginning of column 1 to (3,1)
      {
         Iterator pos = lower.erase( 1UL, lower.begin( 1UL ), lower.find( 3UL, 1UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 7UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(2,0) != 2 || lower(2,2) != 4 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,2) != 7 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 0 4 0 )\n( 5 6 7 8 )\n";
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

      // Erasing the elements from (3,2) to the column end
      {
         Iterator pos = lower.erase( 2UL, lower.find( 3UL, 2UL ), lower.end( 2UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(0,0) != 1 ||
             lower(2,0) != 2 || lower(2,2) != 4 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 2 0 4 0 )\n( 5 6 0 8 )\n";
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

      // Erasing the elements from (0,0) to (3,0)
      {
         Iterator pos = lower.erase( 0UL, lower.find( 0UL, 0UL ), lower.find( 3UL, 0UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(2,2) != 4 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a multi-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 4 0 )\n( 5 6 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = lower.erase( 0UL, lower.find( 3UL, 0UL ), lower.find( 3UL, 0UL ) );

         checkRows    ( lower, 4UL );
         checkColumns ( lower, 4UL );
         checkCapacity( lower, 8UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );
         checkNonZeros( lower, 3UL, 1UL );

         if( lower(2,2) != 4 ||
             lower(3,0) != 5 || lower(3,1) != 6 || lower(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 4 0 )\n( 5 6 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::resize()";

      // Initialization check
      LT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Resizing to 2x2
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );

      // Resizing to 4x4 and preserving the elements
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower.resize( 4UL, true );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );

      // Resizing to 2x2
      lower(2,2) = 4;
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );

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
      test_ = "Column-major LowerMatrix::resize()";

      // Initialization check
      OLT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Resizing to 2x2
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );

      // Resizing to 4x4 and preserving the elements
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower.resize( 4UL, true );

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );

      // Resizing to 2x2
      lower(2,2) = 4;
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );

      // Resizing to 0x0
      lower.resize( 0UL );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::reserve()";

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
      test_ = "Column-major LowerMatrix::reserve()";

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
/*!\brief Test of the \c trim() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::trim()";

      // Initialization check
      LT lower( 3UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 0UL );

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

      // Trimming the matrix
      lower.trim();

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 0UL );
      checkCapacity( lower,  1UL, 0UL );
      checkCapacity( lower,  2UL, 0UL );
   }

   {
      test_ = "Row-major LowerMatrix::trim( size_t )";

      // Initialization check
      LT lower( 3UL, 3UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 0UL );

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

      // Trimming the 0th row
      lower.trim( 0UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL,  0UL );
      checkCapacity( lower,  1UL, 25UL );
      checkCapacity( lower,  2UL, 20UL );

      // Trimming the 1st row
      lower.trim( 1UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL,  0UL );
      checkCapacity( lower,  1UL,  0UL );
      checkCapacity( lower,  2UL, 45UL );

      // Trimming the 2nd row
      lower.trim( 2UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 0UL );
      checkCapacity( lower,  1UL, 0UL );
      checkCapacity( lower,  2UL, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::trim()";

      // Initialization check
      OLT lower( 3UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 0UL );

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

      // Trimming the matrix
      lower.trim();

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 0UL );
      checkCapacity( lower,  1UL, 0UL );
      checkCapacity( lower,  2UL, 0UL );
   }

   {
      test_ = "Column-major LowerMatrix::trim( size_t )";

      // Initialization check
      OLT lower( 3UL, 3UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 0UL );

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

      // Trimming the 0th column
      lower.trim( 0UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL,  0UL );
      checkCapacity( lower,  1UL, 25UL );
      checkCapacity( lower,  2UL, 20UL );

      // Trimming the 1st column
      lower.trim( 1UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL,  0UL );
      checkCapacity( lower,  1UL,  0UL );
      checkCapacity( lower,  2UL, 45UL );

      // Trimming the 2nd column
      lower.trim( 2UL );

      checkRows    ( lower,  3UL );
      checkColumns ( lower,  3UL );
      checkCapacity( lower, 45UL );
      checkCapacity( lower,  0UL, 0UL );
      checkCapacity( lower,  1UL, 0UL );
      checkCapacity( lower,  2UL, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the LowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix swap";

      LT lower1( 2UL );
      lower1(0,0) = 1;
      lower1(1,0) = 2;
      lower1(1,1) = 3;

      LT lower2( 2UL );
      lower2(0,0) = 4;
      lower2(1,0) = 5;
      lower2(1,1) = 0;

      swap( lower1, lower2 );

      checkRows    ( lower1, 2UL );
      checkColumns ( lower1, 2UL );
      checkCapacity( lower1, 2UL );
      checkNonZeros( lower1, 2UL );
      checkNonZeros( lower1, 0UL, 1UL );
      checkNonZeros( lower1, 1UL, 1UL );

      if( lower1(0,0) != 4 || lower1(0,1) != 0 || lower1(1,0) != 5 || lower1(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower1 << "\n"
             << "   Expected result:\n( 4 0 )\n( 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( lower2, 2UL );
      checkColumns ( lower2, 2UL );
      checkCapacity( lower2, 3UL );
      checkNonZeros( lower2, 3UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );

      if( lower2(0,0) != 1 || lower2(0,1) != 0 || lower2(1,0) != 2 || lower2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix swap";

      OLT lower1( 2UL );
      lower1(0,0) = 1;
      lower1(1,0) = 2;
      lower1(1,1) = 3;

      OLT lower2( 2UL );
      lower2(0,0) = 4;
      lower2(1,0) = 5;
      lower2(1,1) = 0;

      swap( lower1, lower2 );

      checkRows    ( lower1, 2UL );
      checkColumns ( lower1, 2UL );
      checkCapacity( lower1, 2UL );
      checkNonZeros( lower1, 2UL );
      checkNonZeros( lower1, 0UL, 2UL );
      checkNonZeros( lower1, 1UL, 0UL );

      if( lower1(0,0) != 4 || lower1(0,1) != 0 || lower1(1,0) != 5 || lower1(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower1 << "\n"
             << "   Expected result:\n( 4 0 )\n( 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( lower2, 2UL );
      checkColumns ( lower2, 2UL );
      checkCapacity( lower2, 3UL );
      checkNonZeros( lower2, 3UL );
      checkNonZeros( lower2, 0UL, 2UL );
      checkNonZeros( lower2, 1UL, 1UL );

      if( lower2(0,0) != 1 || lower2(0,1) != 0 || lower2(1,0) != 2 || lower2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::find()";

      typedef LT::ConstIterator  ConstIterator;

      // Initialization check
      LT lower( 8UL, 3UL );
      lower(2,1) = 1;
      lower(3,2) = 2;
      lower(6,5) = 3;

      checkRows    ( lower, 8UL );
      checkColumns ( lower, 8UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 1UL );
      checkNonZeros( lower, 4UL, 0UL );
      checkNonZeros( lower, 5UL, 0UL );
      checkNonZeros( lower, 6UL, 1UL );
      checkNonZeros( lower, 7UL, 0UL );

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
         else if( pos->index() != 1 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( lower.find( 3UL, 2UL ) );

         if( pos == lower.end( 3UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (3,2)\n"
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

      // Searching for the third element
      {
         ConstIterator pos( lower.find( 6UL, 5UL ) );

         if( pos == lower.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,5)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 5 || pos->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 5\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 3\n"
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
      test_ = "Column-major LowerMatrix::find()";

      typedef OLT::ConstIterator  ConstIterator;

      // Initialization check
      OLT lower( 8UL, 3UL );
      lower(2,1) = 1;
      lower(3,2) = 2;
      lower(6,5) = 3;

      checkRows    ( lower, 8UL );
      checkColumns ( lower, 8UL );
      checkCapacity( lower, 3UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );
      checkNonZeros( lower, 3UL, 0UL );
      checkNonZeros( lower, 4UL, 0UL );
      checkNonZeros( lower, 5UL, 1UL );
      checkNonZeros( lower, 6UL, 0UL );
      checkNonZeros( lower, 7UL, 0UL );

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

      // Searching for the second element
      {
         ConstIterator pos( lower.find( 3UL, 2UL ) );

         if( pos == lower.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (3,2)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 3 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 3\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( lower.find( 6UL, 5UL ) );

         if( pos == lower.end( 5UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,5)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 6 || pos->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 6\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 3\n"
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
/*!\brief Test of the \c lowerBound() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::lowerBound()";

      typedef LT::ConstIterator  ConstIterator;

      // Initialization check
      LT lower( 6UL, 2UL );
      lower(4,1) = 1;
      lower(4,3) = 2;

      checkRows    ( lower, 6UL );
      checkColumns ( lower, 6UL );
      checkCapacity( lower, 2UL );
      checkNonZeros( lower, 2UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );
      checkNonZeros( lower, 4UL, 2UL );
      checkNonZeros( lower, 5UL, 0UL );

      // Determining the lower bound for position (4,0)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 0UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,0)\n"
                << "   Current matrix:\n" << lower << "\n";
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
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

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
         else if( pos->index() != 1 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
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
         else if( pos->index() != 3 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 3\n"
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
         else if( pos->index() != 3 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 3\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,4)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 4UL ) );

         if( pos != lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,4)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::lowerBound()";

      typedef OLT::ConstIterator  ConstIterator;

      // Initialization check
      OLT lower( 6UL, 2UL );
      lower(2,1) = 1;
      lower(4,1) = 2;

      checkRows    ( lower, 6UL );
      checkColumns ( lower, 6UL );
      checkCapacity( lower, 2UL );
      checkNonZeros( lower, 2UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );
      checkNonZeros( lower, 4UL, 0UL );
      checkNonZeros( lower, 5UL, 0UL );

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( lower.lowerBound( 1UL, 1UL ) );

         if( pos == lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
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

      // Determining the lower bound for position (2,1)
      {
         ConstIterator pos( lower.lowerBound( 2UL, 1UL ) );

         if( pos == lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
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

      // Determining the lower bound for position (3,1)
      {
         ConstIterator pos( lower.lowerBound( 3UL, 1UL ) );

         if( pos == lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,1)\n"
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

      // Determining the lower bound for position (4,1)
      {
         ConstIterator pos( lower.lowerBound( 4UL, 1UL ) );

         if( pos == lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
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

      // Determining the lower bound for position (5,1)
      {
         ConstIterator pos( lower.lowerBound( 5UL, 1UL ) );

         if( pos != lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,1)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::upperBound()";

      typedef LT::ConstIterator  ConstIterator;

      // Initialization check
      LT lower( 6UL, 2UL );
      lower(4,1) = 1;
      lower(4,3) = 2;

      checkRows    ( lower, 6UL );
      checkColumns ( lower, 6UL );
      checkCapacity( lower, 2UL );
      checkNonZeros( lower, 2UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );
      checkNonZeros( lower, 4UL, 2UL );
      checkNonZeros( lower, 5UL, 0UL );

      // Determining the upper bound for position (4,0)
      {
         ConstIterator pos( lower.upperBound( 4UL, 0UL ) );

         if( pos == lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,0)\n"
                << "   Current matrix:\n" << lower << "\n";
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
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

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
         else if( pos->index() != 3 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 3\n"
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
         else if( pos->index() != 3 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 3\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,3)
      {
         ConstIterator pos( lower.upperBound( 4UL, 3UL ) );

         if( pos != lower.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,3)\n"
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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::lowerBound()";

      typedef OLT::ConstIterator  ConstIterator;

      // Initialization check
      OLT lower( 6UL, 2UL );
      lower(2,1) = 1;
      lower(4,1) = 2;

      checkRows    ( lower, 6UL );
      checkColumns ( lower, 6UL );
      checkCapacity( lower, 2UL );
      checkNonZeros( lower, 2UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 0UL );
      checkNonZeros( lower, 3UL, 0UL );
      checkNonZeros( lower, 4UL, 0UL );
      checkNonZeros( lower, 5UL, 0UL );

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( lower.upperBound( 1UL, 1UL ) );

         if( pos == lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
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

      // Determining the upper bound for position (2,1)
      {
         ConstIterator pos( lower.upperBound( 2UL, 1UL ) );

         if( pos == lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
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

      // Determining the upper bound for position (3,1)
      {
         ConstIterator pos( lower.upperBound( 3UL, 1UL ) );

         if( pos == lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,1)\n"
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

      // Determining the upper bound for position (4,1)
      {
         ConstIterator pos( lower.upperBound( 4UL, 1UL ) );

         if( pos != lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (5,1)
      {
         ConstIterator pos( lower.upperBound( 5UL, 1UL ) );

         if( pos != lower.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,1)\n"
                << "   Current matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the LowerMatrix
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
/*!\brief Test of the assignment to submatrices of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to submatrices of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      typedef blaze::SparseSubmatrix<LT>  SMT;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      SMT sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      if( sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(1,1) << "\n"
             << "   Expected result: 3\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,0) = -5;

      if( sm(0,0) !=  2 || sm(0,1) != 0 ||
          sm(1,0) != -5 || sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  2  0 )\n( -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4  2  0 )\n(  7 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 0 0 )\n(  7 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      typedef blaze::SparseSubmatrix<OLT>  SMT;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      SMT sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      if( sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(1,1) << "\n"
             << "   Expected result: 3\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,0) = -5;

      if( sm(0,0) !=  2 || sm(0,1) != 0 ||
          sm(1,0) != -5 || sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  2  0 )\n( -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4  2  0 )\n(  7 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 0 0 )\n(  7 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the assignment to rows of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to rows of the LowerMatrix specialization.
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      RT row1 = row( lower, 1UL );

      if( row1[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: 2\n";
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

      row1[1] = -5;

      if( row1[0] != -4 || row1[1] != -5 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -4 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -5 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4 -5  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != 0 || row1[1] != 0 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 7 0 3 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      RT row1 = row( lower, 1UL );

      if( row1[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: 2\n";
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

      row1[1] = -5;

      if( row1[0] != -4 || row1[1] != -5 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -4 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -5 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4 -5  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != 0 || row1[1] != 0 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the assignment to columns of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to columns of the LowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      CT col1 = column( lower, 1UL );

      if( col1[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it = col1.begin();

      if( it == col1.end() || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      col1[1] = -5;

      if( col1[0] != 0 || col1[1] != -5 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -5 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4 -5  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col1 );

      if( col1[0] != 0 || col1[1] != 0 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 0 0 )\n(  7 0 3 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      CT col1 = column( lower, 1UL );

      if( col1[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it = col1.begin();

      if( it == col1.end() || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      col1[1] = -5;

      if( col1[0] != 0 || col1[1] != -5 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -5 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4 -5  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col1 );

      if( col1[0] != 0 || col1[1] != 0 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 0 0 )\n(  7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace lowermatrix

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
   std::cout << "   Running LowerMatrix sparse test..." << std::endl;

   try
   {
      RUN_LOWERMATRIX_SPARSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during LowerMatrix sparse test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
