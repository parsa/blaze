//=================================================================================================
/*!
//  \file src/mathtest/diagonalmatrix/DenseTest.cpp
//  \brief Source file for the DiagonalMatrix dense test
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
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DenseColumn.h>
#include <blaze/math/DenseRow.h>
#include <blaze/math/DenseSubmatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/diagonalmatrix/DenseTest.h>


namespace blazetest {

namespace mathtest {

namespace diagonalmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DiagonalMatrix dense test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseTest::DenseTest()
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
   testResize();
   testExtend();
   testReserve();
   testSwap();
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
/*!\brief Test of the DiagonalMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the DiagonalMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Row-major DiagonalMatrix default constructor (StaticMatrix)";

      const blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > diag;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Row-major DiagonalMatrix default constructor (HybridMatrix)";

      const blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Row-major DiagonalMatrix default constructor (DynamicMatrix)";

      const DT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   // Size constructor (HybridMatrix)
   {
      test_ = "Row-major DiagonalMatrix size constructor (HybridMatrix)";

      const blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > diag( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 0UL );
   }

   // Size constructor (DynamicMatrix)
   {
      test_ = "Row-major DiagonalMatrix size constructor (DynamicMatrix)";

      const DT diag( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Row-major DiagonalMatrix copy constructor (0x0)";

      const DT diag1;
      const DT diag2( diag1 );

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Row-major DiagonalMatrix copy constructor (3x3)";

      DT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      const DT diag2( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Row-major DiagonalMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      const DT diag( mat );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Conversion constructor (diagonal)
   {
      test_ = "Row-major DiagonalMatrix conversion constructor (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      const DT diag( mat );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (lower)
   {
      test_ = "Row-major DiagonalMatrix conversion constructor (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         const DT diag( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-diagonal DiagonalMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (upper)
   {
      test_ = "Row-major DiagonalMatrix conversion constructor (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         const DT diag( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-diagonal DiagonalMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (DiagonalMatrix)
   {
      test_ = "Row-major DiagonalMatrix conversion constructor (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      const DT diag2( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Column-major DiagonalMatrix default constructor (StaticMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > diag;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Column-major DiagonalMatrix default constructor (HybridMatrix)";

      blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Column-major DiagonalMatrix default constructor (DynamicMatrix)";

      ODT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   // Size constructor (HybridMatrix)
   {
      test_ = "Column-major DiagonalMatrix size constructor (HybridMatrix)";

      blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > diag( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkNonZeros( diag, 0UL );
   }

   // Size constructor (DynamicMatrix)
   {
      test_ = "Column-major DiagonalMatrix size constructor (DynamicMatrix)";

      ODT diag( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Column-major DiagonalMatrix copy constructor (0x0)";

      const ODT diag1;
      const ODT diag2( diag1 );

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Column-major DiagonalMatrix copy constructor (3x3)";

      ODT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      const ODT diag2( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Column-major DiagonalMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      const ODT diag( mat );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Conversion constructor (diagonal)
   {
      test_ = "Column-major DiagonalMatrix conversion constructor (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      const ODT diag( mat );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (lower)
   {
      test_ = "Column-major DiagonalMatrix conversion constructor (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         const ODT diag( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-diagonal DiagonalMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (upper)
   {
      test_ = "Column-major DiagonalMatrix conversion constructor (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         const ODT diag( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-diagonal DiagonalMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (DiagonalMatrix)
   {
      test_ = "Column-major DiagonalMatrix conversion constructor (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      const ODT diag2( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the DiagonalMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   // Homogeneous assignment (3x3)
   {
      test_ = "Row-major DiagonalMatrix homogeneous assignment (3x3)";

      DT diag( 3UL );
      diag = 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 2 0 )\n( 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Row-major DiagonalMatrix copy assignment (0x0)";

      DT diag1, diag2;

      diag2 = diag1;

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Row-major DiagonalMatrix copy assignment (3x3)";

      DT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major DiagonalMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      DT diag;
      diag = mat;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Row-major/row-major dense matrix assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix assignment (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      DT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix assignment (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      DT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix assignment (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix assignment (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major DiagonalMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      DT diag;
      diag = mat;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   // Homogeneous assignment (3x3)
   {
      test_ = "Column-major DiagonalMatrix homogeneous assignment (3x3)";

      ODT diag( 3UL );
      diag = 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 2 0 )\n( 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Column-major DiagonalMatrix copy assignment (0x0)";

      ODT diag1, diag2;

      diag2 = diag1;

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Column-major DiagonalMatrix copy assignment (3x3)";

      ODT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major DiagonalMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      ODT diag;
      diag = mat;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Column-major/row-major dense matrix assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix assignment (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      ODT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix assignment (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      ODT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix assignment (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix assignment (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major DiagonalMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      ODT diag;
      diag = mat;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   // Row-major/row-major dense matrix addition assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix addition assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) = -2;
      mat(2,2) =  2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix addition assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) = -2;
      mat(2,2) =  2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix addition assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix addition assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix addition assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix addition assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix addition assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix addition assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix addition assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) = -2;
      mat(2,2) =  2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix addition assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) = -2;
      mat(2,2) =  2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix addition assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix addition assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix addition assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix addition assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix addition assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix addition assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix addition assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix addition assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   // Column-major/row-major dense matrix addition assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix addition assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) = -2;
      mat(2,2) =  2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix addition assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) = -2;
      mat(2,2) =  2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix addition assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix addition assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix addition assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix addition assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix addition assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix addition assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix addition assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) = -2;
      mat(2,2) =  2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix addition assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) = -2;
      mat(2,2) =  2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix addition assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix addition assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix addition assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix addition assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix addition assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix addition assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix addition assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix addition assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major dense matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) =  2;
      mat(2,2) = -2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) =  2;
      mat(2,2) = -2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix subtraction assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix subtraction assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix subtraction assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) =  2;
      mat(2,2) = -2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) =  2;
      mat(2,2) = -2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix subtraction assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix subtraction assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix subtraction assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix subtraction assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix subtraction assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix subtraction assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major dense matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) =  2;
      mat(2,2) = -2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) =  2;
      mat(2,2) = -2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix subtraction assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix subtraction assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix subtraction assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix subtraction assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix subtraction assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix subtraction assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) =  2;
      mat(2,2) = -2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) =  2;
      mat(2,2) = -2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix subtraction assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix subtraction assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix subtraction assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix subtraction assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix subtraction assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix subtraction assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major dense matrix multiplication assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix multiplication assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix multiplication assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix multiplication assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix multiplication assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix multiplication assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix multiplication assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix multiplication assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix multiplication assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix multiplication assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix multiplication assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix multiplication assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix multiplication assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix multiplication assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix multiplication assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major dense matrix multiplication assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix multiplication assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix multiplication assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix multiplication assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (upper)
   {
      test_ = "Column";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix multiplication assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix multiplication assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix multiplication assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix multiplication assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix multiplication assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix multiplication assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix multiplication assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix multiplication assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix multiplication assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix multiplication assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all DiagonalMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      DT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag *= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      DT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag = diag * 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      DT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag = 2 * diag;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      DT diag( 3UL );
      diag(1,1) =  4;
      diag(2,2) = -6;

      diag /= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      DT diag( 3UL );
      diag(1,1) =  4;
      diag(2,2) = -6;

      diag = diag / 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major DiagonalMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::scale()";

      // Initialization check
      DT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      diag.scale( 2 );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      diag.scale( 0.5 );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major DiagonalMatrix::scale() (complex)";

      using blaze::complex;

      blaze::DiagonalMatrix< blaze::DynamicMatrix<complex<float>,blaze::rowMajor> > diag( 2UL );
      diag(0,0) = complex<float>( 1.0F, 0.0F );
      diag(1,1) = complex<float>( 2.0F, 0.0F );

      diag.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );

      if( diag(0,0) != complex<float>( 3.0F, 0.0F ) || diag(0,1) != complex<float>( 0.0F, 0.0F ) ||
          diag(1,0) != complex<float>( 0.0F, 0.0F ) || diag(1,1) != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( (3,0) (0,0)\n(0,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      ODT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag *= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      ODT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag = diag * 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      ODT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag = 2 * diag;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      ODT diag( 3UL );
      diag(1,1) =  4;
      diag(2,2) = -6;

      diag /= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      ODT diag( 3UL );
      diag(1,1) =  4;
      diag(2,2) = -6;

      diag = diag / 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major DiagonalMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::scale()";

      // Initialization check
      ODT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      diag.scale( 2 );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      diag.scale( 0.5 );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major DiagonalMatrix::scale() (complex)";

      using blaze::complex;

      blaze::DiagonalMatrix< blaze::DynamicMatrix<complex<float>,blaze::columnMajor> > diag( 2UL );
      diag(0,0) = complex<float>( 1.0F, 0.0F );
      diag(1,1) = complex<float>( 2.0F, 0.0F );

      diag.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );

      if( diag(0,0) != complex<float>( 3.0F, 0.0F ) || diag(0,1) != complex<float>( 0.0F, 0.0F ) ||
          diag(1,0) != complex<float>( 0.0F, 0.0F ) || diag(1,1) != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( (3,0) (0,0)\n(0,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the DiagonalMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void DenseTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::operator()";

      DT diag( 3UL );

      // Writing the element (1,1)
      diag(1,1) = 1;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the element (2,2)
      diag(2,2) = 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Adding to the element (0,0)
      diag(0,0) += 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 3 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3 0 0 )\n( 0 1 0 )\n( 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtracting from the element (1,1)
      diag(1,1) -= 4;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0  2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplying the element (2,2)
      diag(2,2) *= -3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Dividing the element (2,2)
      diag(2,2) /= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to write the element (2,1)
      try {
         diag(2,1) = 5;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to lower matrix element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to write the element (1,2)
      try {
         diag(1,2) = 5;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to upper matrix element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::operator()";

      ODT diag( 3UL );

      // Writing the element (1,1)
      diag(1,1) = 1;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the element (2,2)
      diag(2,2) = 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Adding to the element (0,0)
      diag(0,0) += 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 3 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3 0 0 )\n( 0 1 0 )\n( 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtracting from the element (1,1)
      diag(1,1) -= 4;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0  2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplying the element (2,2)
      diag(2,2) *= -3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Dividing the element (2,2)
      diag(2,2) /= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to write the element (2,1)
      try {
         diag(2,1) = 5;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to lower matrix element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to write the element (1,2)
      try {
         diag(1,2) = 5;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to upper matrix element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      typedef DT::Iterator       Iterator;
      typedef DT::ConstIterator  ConstIterator;

      DT diag( 3UL );
      diag(0,0) =  1;
      diag(1,1) = -2;
      diag(2,2) =  3;

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Row-major Iterator/ConstIterator conversion";

         ConstIterator it( begin( diag, 1UL ) );

         if( it == end( diag, 1UL ) || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator
      {
         test_ = "Row-major Iterator subtraction";

         const size_t number( end( diag, 0UL ) - begin( diag, 0UL ) );

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

         const size_t number( cend( diag, 1UL ) - cbegin( diag, 1UL ) );

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

         ConstIterator it ( cbegin( diag, 2UL ) );
         ConstIterator end( cend( diag, 2UL ) );

         if( it == end || *it != 0 ) {
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

         if( it == end || *it != 0 ) {
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

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 3 ) {
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

         if( it == end || *it != 3 ) {
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

      // Testing assignment to diagonal elements via Iterator
      {
         test_ = "Row-major assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it = 4;

         if( diag(0,0) != 4 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4  0  0 )\n( 0 -2  0 )\n( 0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Row-major assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 1UL );
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing assignment to upper elements via Iterator
      {
         test_ = "Row-major assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 0UL ) + 1UL;
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to diagonal elements via Iterator
      {
         test_ = "Row-major addition assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 1UL ) + 1UL;
         *it += 3;

         if( diag(0,0) != 4 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4 0 0 )\n( 0 1 0 )\n( 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Row-major addition assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 2UL );
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to upper elements via Iterator
      {
         test_ = "Row-major addition assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 0UL ) + 2UL;
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to diagonal elements via Iterator
      {
         test_ = "Row-major subtraction assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 2UL ) + 2UL;
         *it -= 4;

         if( diag(0,0) != 4 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Row-major subtraction assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 2UL ) + 1UL;
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to upper elements via Iterator
      {
         test_ = "Row-major subtraction assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 1UL ) + 2UL;
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to diagonal elements via Iterator
      {
         test_ = "Row-major multiplication assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it *= 2;

         if( diag(0,0) != 8 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 8  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Row-major multiplication assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 1UL );
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to upper elements via Iterator
      {
         test_ = "Row-major multiplication assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 0UL ) + 1UL;
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to diagonal elements via Iterator
      {
         test_ = "Row-major division assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it /= 4;

         if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 2  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Row-major division assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 2UL );
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to upper elements via Iterator
      {
         test_ = "Row-major division assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 0UL ) + 2UL;
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      typedef ODT::Iterator       Iterator;
      typedef ODT::ConstIterator  ConstIterator;

      ODT diag( 3UL );
      diag(0,0) =  1;
      diag(1,1) = -2;
      diag(2,2) =  3;

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Row-major Iterator/ConstIterator conversion";

         ConstIterator it( begin( diag, 1UL ) );

         if( it == end( diag, 1UL ) || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator
      {
         test_ = "Row-major Iterator subtraction";

         const size_t number( end( diag, 0UL ) - begin( diag, 0UL ) );

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

         const size_t number( cend( diag, 1UL ) - cbegin( diag, 1UL ) );

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

         ConstIterator it ( cbegin( diag, 2UL ) );
         ConstIterator end( cend( diag, 2UL ) );

         if( it == end || *it != 0 ) {
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

         if( it == end || *it != 0 ) {
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

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 3 ) {
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

         if( it == end || *it != 3 ) {
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

      // Testing assignment to diagonal elements via Iterator
      {
         test_ = "Row-major assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it = 4;

         if( diag(0,0) != 4 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4  0  0 )\n( 0 -2  0 )\n( 0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Row-major assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 1UL );
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing assignment to upper elements via Iterator
      {
         test_ = "Row-major assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 0UL ) + 1UL;
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to diagonal elements via Iterator
      {
         test_ = "Row-major addition assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 1UL ) + 1UL;
         *it += 3;

         if( diag(0,0) != 4 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4 0 0 )\n( 0 1 0 )\n( 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Row-major addition assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 2UL );
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to upper elements via Iterator
      {
         test_ = "Row-major addition assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 0UL ) + 2UL;
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to diagonal elements via Iterator
      {
         test_ = "Row-major subtraction assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 2UL ) + 2UL;
         *it -= 4;

         if( diag(0,0) != 4 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Row-major subtraction assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 2UL ) + 1UL;
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to upper elements via Iterator
      {
         test_ = "Row-major subtraction assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 1UL ) + 2UL;
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to diagonal elements via Iterator
      {
         test_ = "Row-major multiplication assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it *= 2;

         if( diag(0,0) != 8 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 8  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Row-major multiplication assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 1UL );
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to upper elements via Iterator
      {
         test_ = "Row-major multiplication assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 0UL ) + 1UL;
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to diagonal elements via Iterator
      {
         test_ = "Row-major division assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it /= 4;

         if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 2  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Row-major division assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( diag, 2UL );
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to upper elements via Iterator
      {
         test_ = "Row-major division assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( diag, 0UL ) + 2UL;
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::nonZeros()";

      // Empty matrix
      {
         DT diag( 3UL );

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 9UL );
         checkNonZeros( diag, 0UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         DT diag( 3UL );
         diag(0,0) =  1;
         diag(1,1) = -2;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 9UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         DT diag( 3UL );
         diag(0,0) = -1;
         diag(1,1) =  2;
         diag(2,2) =  3;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 9UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != -1 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) !=  0 || diag(1,1) != 2 || diag(1,2) != 0 ||
             diag(2,0) !=  0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( -1  0  0 )\n(  0  2  0 )\n(  0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::nonZeros()";

      // Empty matrix
      {
         ODT diag( 3UL );

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 9UL );
         checkNonZeros( diag, 0UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         ODT diag( 3UL );
         diag(0,0) =  1;
         diag(1,1) = -2;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 9UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         ODT diag( 3UL );
         diag(0,0) = -1;
         diag(1,1) =  2;
         diag(2,2) =  3;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 9UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != -1 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) !=  0 || diag(1,1) != 2 || diag(1,2) != 0 ||
             diag(2,0) !=  0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( -1  0  0 )\n(  0  2  0 )\n(  0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::reset()";

      // Initialization check
      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a diagonal element
      reset( diag(1,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( diag(1,0) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( diag(0,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 2
      reset( diag, 2UL );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( diag );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::reset()";

      // Initialization check
      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a diagonal element
      reset( diag(1,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( diag(1,0) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( diag(0,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 2
      reset( diag, 2UL );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( diag );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::clear()";

      // Initialization check
      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( diag(1,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( diag(1,0) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( diag(0,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( diag );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::clear()";

      // Initialization check
      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( diag(1,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( diag(1,0) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( diag(0,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( diag );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::resize()";

      // Initialization check
      DT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Resizing to 2x2
      diag.resize( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );

      if( diag(0,1) != 0 || diag(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( x 0 )\n( 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag.resize( 4UL, true );

      checkRows    ( diag,  4UL );
      checkColumns ( diag,  4UL );
      checkCapacity( diag, 16UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 x 0 )\n( 0 0 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      diag(2,2) = 3;
      diag.resize( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      diag.resize( 0UL );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::resize()";

      // Initialization check
      ODT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Resizing to 2x2
      diag.resize( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );

      if( diag(0,1) != 0 || diag(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( x 0 )\n( 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag.resize( 4UL, true );

      checkRows    ( diag,  4UL );
      checkColumns ( diag,  4UL );
      checkCapacity( diag, 16UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 x 0 )\n( 0 0 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      diag(2,2) = 3;
      diag.resize( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      diag.resize( 0UL );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c extend() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::extend()";

      // Initialization check
      DT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Extending the size of the matrix to 2x2
      diag.extend( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );

      if( diag(0,1) != 0 || diag(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( x 0 )\n( 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag.extend( 2UL, true );

      checkRows    ( diag,  4UL );
      checkColumns ( diag,  4UL );
      checkCapacity( diag, 16UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 x 0 )\n( 0 0 0 x )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::extend()";

      // Initialization check
      ODT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Extending the size of the matrix to 2x2
      diag.extend( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );

      if( diag(0,1) != 0 || diag(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( x 0 )\n( 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag.extend( 2UL, true );

      checkRows    ( diag,  4UL );
      checkColumns ( diag,  4UL );
      checkCapacity( diag, 16UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 x 0 )\n( 0 0 0 x )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::reserve()";

      // Initialization check
      DT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Increasing the capacity of the matrix
      diag.reserve( 10UL );

      checkRows    ( diag,  0UL );
      checkColumns ( diag,  0UL );
      checkCapacity( diag, 10UL );
      checkNonZeros( diag,  0UL );

      // Further increasing the capacity of the matrix
      diag.reserve( 20UL );

      checkRows    ( diag,  0UL );
      checkColumns ( diag,  0UL );
      checkCapacity( diag, 20UL );
      checkNonZeros( diag,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::reserve()";

      // Initialization check
      ODT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Increasing the capacity of the matrix
      diag.reserve( 10UL );

      checkRows    ( diag,  0UL );
      checkColumns ( diag,  0UL );
      checkCapacity( diag, 10UL );
      checkNonZeros( diag,  0UL );

      // Further increasing the capacity of the matrix
      diag.reserve( 20UL );

      checkRows    ( diag,  0UL );
      checkColumns ( diag,  0UL );
      checkCapacity( diag, 20UL );
      checkNonZeros( diag,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the DiagonalMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix swap";

      DT diag1( 2UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 3;
      diag2(1,1) = 4;
      diag2(2,2) = 5;

      swap( diag1, diag2 );

      checkRows    ( diag1, 3UL );
      checkColumns ( diag1, 3UL );
      checkCapacity( diag1, 9UL );
      checkNonZeros( diag1, 3UL );
      checkNonZeros( diag1, 0UL, 1UL );
      checkNonZeros( diag1, 1UL, 1UL );
      checkNonZeros( diag1, 2UL, 1UL );

      if( diag1(0,0) != 3 || diag1(0,1) != 0 || diag1(0,2) != 0 ||
          diag1(1,0) != 0 || diag1(1,1) != 4 || diag1(1,2) != 0 ||
          diag1(2,0) != 0 || diag1(2,1) != 0 || diag1(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag1 << "\n"
             << "   Expected result:\n( 3 0 0 )\n( 0 4 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( diag2, 2UL );
      checkColumns ( diag2, 2UL );
      checkCapacity( diag2, 4UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(1,0) != 0 || diag2(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix swap";

      ODT diag1( 2UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;

      ODT diag2( 3UL );
      diag2(0,0) = 3;
      diag2(1,1) = 4;
      diag2(2,2) = 5;

      swap( diag1, diag2 );

      checkRows    ( diag1, 3UL );
      checkColumns ( diag1, 3UL );
      checkCapacity( diag1, 9UL );
      checkNonZeros( diag1, 3UL );
      checkNonZeros( diag1, 0UL, 1UL );
      checkNonZeros( diag1, 1UL, 1UL );
      checkNonZeros( diag1, 2UL, 1UL );

      if( diag1(0,0) != 3 || diag1(0,1) != 0 || diag1(0,2) != 0 ||
          diag1(1,0) != 0 || diag1(1,1) != 4 || diag1(1,2) != 0 ||
          diag1(2,0) != 0 || diag1(2,1) != 0 || diag1(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag1 << "\n"
             << "   Expected result:\n( 3 0 0 )\n( 0 4 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( diag2, 2UL );
      checkColumns ( diag2, 2UL );
      checkCapacity( diag2, 4UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(1,0) != 0 || diag2(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         DT diag;

         if( isDefault( diag ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         DT diag( 3UL );

         if( isDefault( diag(1,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << diag(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( diag ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         DT diag( 3UL );
         diag(1,1) = 1;

         if( isDefault( diag(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << diag(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( diag ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
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
         ODT diag;

         if( isDefault( diag ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         ODT diag( 3UL );

         if( isDefault( diag(1,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << diag(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( diag ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         ODT diag( 3UL );
         diag(1,1) = 1;

         if( isDefault( diag(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << diag(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( diag ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c submatrix() function with the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c submatrix() function with the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

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

      if( it == sm.end(0UL) || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      sm(0,0) = -5;

      if( sm(0,0) != -5 || sm(0,1) != 0 ||
          sm(1,0) !=  0 || sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -5  0 )\n(  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 1)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = 12;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 12 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 12 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  0 )\n(  0 12 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 2)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = 12;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 12  0  0 )\n( 0  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  3  0 )      ( 0  0  3  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 3)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 2UL, 2UL, 2UL );
      sm = 12;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0  2  0  0 )\n"
                                     "( 0  0  3  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  3  0 )      ( 0  0  3  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 4)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 2UL, 0UL, 2UL, 2UL );
      sm = 12;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0  2  0  0 )\n"
                                     "( 0  0  3  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (row-major dense matrix assignment test 1)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,1) = 11;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (row-major dense matrix assignment test 2)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 18;
      mat(1,2) = 11;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major submatrix() function (row-major dense matrix assignment test 3)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major submatrix() function (row-major dense matrix assignment test 4)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (column-major dense matrix assignment test 1)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,1) = 11;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (column-major dense matrix assignment test 2)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 18;
      mat(1,2) = 11;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (column-major dense matrix assignment test 3)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (column-major dense matrix assignment test 4)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (row-major sparse matrix assignment test 1)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = 18;
      mat(2,1) = 11;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (row-major sparse matrix assignment test 2)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 18;
      mat(1,2) = 11;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (row-major sparse matrix assignment test 3)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;
      mat.insert( 1UL, 0UL, 0 );

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (row-major sparse matrix assignment test 4)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;
      mat.insert( 0UL, 1UL, 0 );

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (column-major sparse matrix assignment test 1)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = 18;
      mat(2,1) = 11;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (column-major sparse matrix assignment test 2)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 18;
      mat(1,2) = 11;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (column-major sparse matrix assignment test 3)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;
      mat.insert( 1UL, 0UL, 0 );

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major submatrix() function (column-major sparse matrix assignment test 4)";

      typedef blaze::DenseSubmatrix<DT>  SMT;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;
      mat.insert( 0UL, 1UL, 0 );

      DT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

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

      if( it == sm.end(0UL) || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      sm(0,0) = -5;

      if( sm(0,0) != -5 || sm(0,1) != 0 ||
          sm(1,0) !=  0 || sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -5  0 )\n(  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 1)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = 12;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 12 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 12 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  0 )\n(  0 12 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 2)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = 12;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 12  0  0 )\n( 0  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  3  0 )      ( 0  0  3  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 3)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 2UL, 2UL, 2UL );
      sm = 12;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0  2  0  0 )\n"
                                     "( 0  0  3  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  3  0 )      ( 0  0  3  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 4)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 2UL, 0UL, 2UL, 2UL );
      sm = 12;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0  2  0  0 )\n"
                                     "( 0  0  3  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (row-major dense matrix assignment test 1)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,1) = 11;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (row-major dense matrix assignment test 2)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 18;
      mat(1,2) = 11;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Column-major submatrix() function (row-major dense matrix assignment test 3)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Column-major submatrix() function (row-major dense matrix assignment test 4)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (column-major dense matrix assignment test 1)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,1) = 11;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (column-major dense matrix assignment test 2)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 18;
      mat(1,2) = 11;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (column-major dense matrix assignment test 3)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (column-major dense matrix assignment test 4)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (row-major sparse matrix assignment test 1)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = 18;
      mat(2,1) = 11;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (row-major sparse matrix assignment test 2)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 18;
      mat(1,2) = 11;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (row-major sparse matrix assignment test 3)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;
      mat.insert( 1UL, 0UL, 0 );

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (row-major sparse matrix assignment test 4)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;
      mat.insert( 0UL, 1UL, 0 );

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (column-major sparse matrix assignment test 1)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = 18;
      mat(2,1) = 11;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (column-major sparse matrix assignment test 2)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 18;
      mat(1,2) = 11;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (column-major sparse matrix assignment test 3)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;
      mat.insert( 1UL, 0UL, 0 );

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major submatrix() function (column-major sparse matrix assignment test 4)";

      typedef blaze::DenseSubmatrix<ODT>  SMT;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;
      mat.insert( 0UL, 1UL, 0 );

      ODT diag( 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c row() function with the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c row() function with the DiagonalMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testRow()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      typedef blaze::DenseRow<DT>  RT;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

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

      if( it == row1.end() || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }

      row1[1] = -5;

      if( row1[0] != 0 || row1[1] != -5 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar assignment
   //=====================================================================================

   {
      test_ = "Row-major row() function (scalar assignment test)";

      typedef blaze::DenseRow<DT>  RT;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );
      row1 = 8;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Row-major row() function (dense vector assignment test 1)";

      typedef blaze::DenseRow<DT>  RT;

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[1] = 8;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );
      row1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 2  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Row-major row() function (dense vector assignment test 2)";

      typedef blaze::DenseRow<DT>  RT;

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;
      vec[1] = 8;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  9 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Row-major row() function (dense vector assignment test 3)";

      typedef blaze::DenseRow<DT>  RT;

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[1] = 8;
      vec[2] = 9;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Row-major row() function (sparse vector assignment test 1)";

      typedef blaze::DenseRow<DT>  RT;

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec.insert( 0UL, 0 );
      vec.insert( 2UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );
      row1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 2 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 2  8  0 )
   // ( 0  0  3 )      ( 7  0  3 )
   {
      test_ = "Row-major row() function (sparse vector assignment test 2)";

      typedef blaze::DenseRow<DT>  RT;

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[0] = 2;
      vec[1] = 8;
      vec.insert( 2UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  9 )
   // ( 0  0  3 )      ( 7  0  3 )
   {
      test_ = "Row-major row() function (sparse vector assignment test 3)";

      typedef blaze::DenseRow<DT>  RT;

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec[2] = 9;
      vec.insert( 0UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      typedef blaze::DenseRow<ODT>  RT;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

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

      if( it == row1.end() || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }

      row1[1] = -5;

      if( row1[0] != 0 || row1[1] != -5 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -4 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar assignment
   //=====================================================================================

   {
      test_ = "Column-major row() function (scalar assignment test)";

      typedef blaze::DenseRow<ODT>  RT;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );
      row1 = 8;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Column-major row() function (dense vector assignment test 1)";

      typedef blaze::DenseRow<ODT>  RT;

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[1] = 8;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );
      row1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 2  8  0 )
   // ( 0  0  3 )      ( 7  0  3 )
   {
      test_ = "Column-major row() function (dense vector assignment test 2)";

      typedef blaze::DenseRow<ODT>  RT;

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;
      vec[1] = 8;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  9 )
   // ( 0  0  3 )      ( 7  0  3 )
   {
      test_ = "Column-major row() function (dense vector assignment test 3)";

      typedef blaze::DenseRow<ODT>  RT;

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[1] = 8;
      vec[2] = 9;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Column-major row() function (sparse vector assignment test 1)";

      typedef blaze::DenseRow<ODT>  RT;

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec.insert( 0UL, 0 );
      vec.insert( 2UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );
      row1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Column-major row() function (sparse vector assignment test 2)";

      typedef blaze::DenseRow<ODT>  RT;

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[0] = 2;
      vec[1] = 8;
      vec.insert( 2UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  9 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Column-major row() function (sparse vector assignment test 3)";

      typedef blaze::DenseRow<ODT>  RT;

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec[2] = 9;
      vec.insert( 0UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c column() function with the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c column() function with the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testColumn()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      typedef blaze::DenseColumn<DT>  CT;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

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

      if( it == col1.end() || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 0\n";
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

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar assignment
   //=====================================================================================

   {
      test_ = "Row-major column() function (scalar assignment test)";

      typedef blaze::DenseColumn<DT>  CT;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );
      col1 = 8;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Row-major column() function (dense vector assignment test 1)";

      typedef blaze::DenseColumn<DT>  CT;

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = 8;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );
      col1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Row-major column() function (dense vector assignment test 2)";

      typedef blaze::DenseColumn<DT>  CT;

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[0] = 9;
      vec[1] = 8;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Row-major column() function (dense vector assignment test 3)";

      typedef blaze::DenseColumn<DT>  CT;

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = 8;
      vec[2] = 2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Row-major column() function (sparse vector assignment test 1)";

      typedef blaze::DenseColumn<DT>  CT;

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec.insert( 0UL, 0 );
      vec.insert( 2UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );
      col1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Row-major column() function (sparse vector assignment test 2)";

      typedef blaze::DenseColumn<DT>  CT;

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = 9;
      vec[1] = 8;
      vec.insert( 2UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Row-major column() function (sparse vector assignment test 3)";

      typedef blaze::DenseColumn<DT>  CT;

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec[2] = 2;
      vec.insert( 0UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      typedef blaze::DenseColumn<ODT>  CT;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

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

      if( it == col1.end() || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 0\n";
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

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar assignment
   //=====================================================================================

   {
      test_ = "Column-major column() function (scalar assignment test)";

      typedef blaze::DenseColumn<ODT>  CT;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );
      col1 = 8;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Column-major column() function (dense vector assignment test 1)";

      typedef blaze::DenseColumn<ODT>  CT;

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = 8;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );
      col1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Column-major column() function (dense vector assignment test 2)";

      typedef blaze::DenseColumn<ODT>  CT;

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[0] = 9;
      vec[1] = 8;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Column-major column() function (dense vector assignment test 3)";

      typedef blaze::DenseColumn<ODT>  CT;

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = 8;
      vec[2] = 2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Column-major column() function (sparse vector assignment test 1)";

      typedef blaze::DenseColumn<ODT>  CT;

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec.insert( 0UL, 0 );
      vec.insert( 2UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );
      col1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Column-major column() function (sparse vector assignment test 2)";

      typedef blaze::DenseColumn<ODT>  CT;

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = 9;
      vec[1] = 8;
      vec.insert( 2UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Column-major column() function (sparse vector assignment test 3)";

      typedef blaze::DenseColumn<ODT>  CT;

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec[2] = 2;
      vec.insert( 0UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************

} // namespace diagonalmatrix

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
   std::cout << "   Running DiagonalMatrix dense test..." << std::endl;

   try
   {
      RUN_DIAGONALMATRIX_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during DiagonalMatrix dense test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
