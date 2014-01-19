//=================================================================================================
/*!
//  \file src/mathtest/sparserow/ClassTest.cpp
//  \brief Source file for the SparseRow class test
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
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Views.h>
#include <blazetest/mathtest/sparserow/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace sparserow {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SparseRow class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
   : mat_ ( 5UL, 4UL )
   , tmat_( 5UL, 4UL )
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
   testReserve();
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
/*!\brief Test of the SparseRow constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the SparseRow class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow constructor";

      initialize();

      // 0th matrix row
      {
         RT row0 = row( mat_, 0UL );

         checkSize    ( row0, 4UL );
         checkNonZeros( row0, 0UL );

         if( row0[0] != 0 || row0[1] != 0 || row0[2] != 0 || row0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 0th sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st matrix row
      {
         RT row1 = row( mat_, 1UL );

         checkSize    ( row1, 4UL );
         checkNonZeros( row1, 1UL );

         if( row1[0] != 0 || row1[1] != 1 || row1[2] != 0 || row1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row1 << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd matrix row
      {
         RT row2 = row( mat_, 2UL );

         checkSize    ( row2, 4UL );
         checkNonZeros( row2, 2UL );

         if( row2[0] != -2 || row2[1] != 0 || row2[2] != -3 || row2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row2 << "\n"
                << "   Expected result:\n( -2 0 -3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd matrix row
      {
         RT row3 = row( mat_, 3UL );

         checkSize    ( row3, 4UL );
         checkNonZeros( row3, 3UL );

         if( row3[0] != 0 || row3[1] != 4 || row3[2] != 5 || row3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th matrix row
      {
         RT row4 = row( mat_, 4UL );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 4UL );

         if( row4[0] != 7 || row4[1] != -8 || row4[2] != 9 || row4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 4th sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 7 -8 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow constructor";

      initialize();

      // 0th matrix row
      {
         TRT row0 = row( tmat_, 0UL );

         checkSize    ( row0, 4UL );
         checkNonZeros( row0, 0UL );

         if( row0[0] != 0 || row0[1] != 0 || row0[2] != 0 || row0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 0th sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st matrix row
      {
         TRT row1 = row( tmat_, 1UL );

         checkSize    ( row1, 4UL );
         checkNonZeros( row1, 1UL );

         if( row1[0] != 0 || row1[1] != 1 || row1[2] != 0 || row1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row1 << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd matrix row
      {
         TRT row2 = row( tmat_, 2UL );

         checkSize    ( row2, 4UL );
         checkNonZeros( row2, 2UL );

         if( row2[0] != -2 || row2[1] != 0 || row2[2] != -3 || row2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row2 << "\n"
                << "   Expected result:\n( -2 0 -3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd matrix row
      {
         TRT row3 = row( tmat_, 3UL );

         checkSize    ( row3, 4UL );
         checkNonZeros( row3, 3UL );

         if( row3[0] != 0 || row3[1] != 4 || row3[2] != 5 || row3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th matrix row
      {
         TRT row4 = row( tmat_, 4UL );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 4UL );

         if( row4[0] != 7 || row4[1] != -8 || row4[2] != 9 || row4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 4th sparse row failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 7 -8 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseRow assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the SparseRow class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseRow copy assignment";

      initialize();

      RT row1 = row( mat_, 1UL );
      row1 = row( mat_, 2UL );

      checkSize    ( row1,  4UL );
      checkNonZeros( row1,  2UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( row1[0] != -2 || row1[1] != 0 || row1[2] != -3 || row1[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -2 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != -2 || mat_(1,1) !=  0 || mat_(1,2) != -3 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector assignment";

      initialize();

      RT row1 = row( mat_, 1UL );

      blaze::DynamicVector<int,blaze::rowVector> vec1( 4UL, 0 );
      vec1[1] = 8;
      vec1[3] = 9;

      row1 = vec1;

      checkSize    ( row1,  4UL );
      checkNonZeros( row1,  2UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 0 || row1[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  8 || mat_(1,2) !=  0 || mat_(1,3) !=  9 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  8  0  9 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector assignment";

      initialize();

      RT row4 = row( mat_, 4UL );

      blaze::CompressedVector<int,blaze::rowVector> vec1( 4UL );
      vec1[3] = 9;

      row4 = vec1;

      checkSize    ( row4, 4UL );
      checkNonZeros( row4, 1UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 7UL );

      if( row4[0] != 0 || row4[1] != 0 || row4[2] != 0 || row4[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  0 || mat_(4,1) != 0 || mat_(4,2) !=  0 || mat_(4,3) !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  0  0  0  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseRow copy assignment";

      initialize();

      TRT row1 = row( tmat_, 1UL );
      row1 = row( tmat_, 2UL );

      checkSize    ( row1 ,  4UL );
      checkNonZeros( row1 ,  2UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( row1[0] != -2 || row1[1] != 0 || row1[2] != -3 || row1[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -2 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != -2 || tmat_(1,1) !=  0 || tmat_(1,2) != -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector assignment";

      initialize();

      TRT row1 = row( tmat_, 1UL );

      blaze::DynamicVector<int,blaze::rowVector> vec1( 4UL, 0 );
      vec1[1] = 8;
      vec1[3] = 9;

      row1 = vec1;

      checkSize    ( row1 ,  4UL );
      checkNonZeros( row1 ,  2UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 0 || row1[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  8 || tmat_(1,2) !=  0 || tmat_(1,3) !=  9 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  8  0  9 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector assignment";

      initialize();

      TRT row4 = row( tmat_, 4UL );

      blaze::CompressedVector<int,blaze::rowVector> vec1( 4UL );
      vec1[3] = 9;

      row4 = vec1;

      checkSize    ( row4 , 4UL );
      checkNonZeros( row4 , 1UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 7UL );

      if( row4[0] != 0 || row4[1] != 0 || row4[2] != 0 || row4[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != 4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  0 || tmat_(4,1) != 0 || tmat_(4,2) !=  0 || tmat_(4,3) !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  0  0  0  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseRow addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the SparseRow class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Row-major SparseRow addition assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseRow addition assignment";

      initialize();

      RT row2 = row( mat_, 2UL );
      row2 += row( mat_, 3UL );

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( row2[0] != -2 || row2[1] != 4 || row2[2] != 2 || row2[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 4 2 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  4 || mat_(2,2) != 2 || mat_(2,3) != -6 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  4  2 -6 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector addition assignment";

      initialize();

      RT row2 = row( mat_, 2UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      row2 += vec;

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  2UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( row2[0] != 0 || row2[1] != -4 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) != -4 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0 -4 -3  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector addition assignment";

      initialize();

      RT row2 = row( mat_, 2UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      row2 += vec;

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  3UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( row2[0] != 0 || row2[1] != -4 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) != -4 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0 -4 -3  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major SparseRow addition assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseRow addition assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );
      row2 += row( tmat_, 3UL );

      checkSize    ( row2 ,  4UL );
      checkNonZeros( row2 ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( row2[0] != -2 || row2[1] != 4 || row2[2] != 2 || row2[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 4 2 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) != 0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  4 || tmat_(2,2) != 2 || tmat_(2,3) != -6 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) != 5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) != 9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  4  2 -6 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector addition assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      row2 += vec;

      checkSize    ( row2 ,  4UL );
      checkNonZeros( row2 ,  2UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( row2[0] != 0 || row2[1] != -4 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 0 || tmat_(2,1) != -4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) != 7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0 -4 -3  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector addition assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      row2 += vec;

      checkSize    ( row2 ,  4UL );
      checkNonZeros( row2 ,  2UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( row2[0] != 0 || row2[1] != -4 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 0 || tmat_(2,1) != -4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) != 7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0 -4 -3  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseRow subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the SparseRow class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Row-major SparseRow subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseRow subtraction assignment";

      initialize();

      RT row2 = row( mat_, 2UL );
      row2 -= row( mat_, 3UL );

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( row2[0] != -2 || row2[1] != -4 || row2[2] != -8 || row2[3] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 -4 -8 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) != -4 || mat_(2,2) != -8 || mat_(2,3) !=  6 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2 -4 -8  6 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector subtraction assignment";

      initialize();

      RT row2 = row( mat_, 2UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      row2 -= vec;

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  3UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( row2[0] != -4 || row2[1] != 4 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=  4 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  4 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector subtraction assignment";

      initialize();

      RT row2 = row( mat_, 2UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      row2 -= vec;

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  3UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( row2[0] != -4 || row2[1] != 4 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=  4 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  4 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major SparseRow subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseRow subtraction assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );
      row2 -= row( tmat_, 3UL );

      checkSize    ( row2 ,  4UL );
      checkNonZeros( row2 ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( row2[0] != -2 || row2[1] != -4 || row2[2] != -8 || row2[3] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 -4 -8 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) != -4 || tmat_(2,2) != -8 || tmat_(2,3) !=  6 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2 -4 -8  6 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector subtraction assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      row2 -= vec;

      checkSize    ( row2 ,  4UL );
      checkNonZeros( row2 ,  3UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( row2[0] != -4 || row2[1] != 4 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -4 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  4 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector subtraction assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      row2 -= vec;

      checkSize    ( row2 ,  4UL );
      checkNonZeros( row2 ,  3UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( row2[0] != -4 || row2[1] != 4 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -4 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  4 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseRow multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the SparseRow class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Row-major SparseRow multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseRow multiplication assignment";

      initialize();

      RT row2 = row( mat_, 2UL );
      row2 *= row( mat_, 3UL );

      checkSize    ( row2, 4UL );
      checkNonZeros( row2, 1UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( row2[0] != 0 || row2[1] != 0 || row2[2] != -15 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( 0 0 -15 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) !=  1 || mat_(1,2) !=   0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) !=  0 || mat_(2,2) != -15 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  4 || mat_(3,2) !=   5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0   0  0 )\n"
                                     "(  0  1   0  0 )\n"
                                     "(  0  0 -15  0 )\n"
                                     "(  0  4   5 -6 )\n"
                                     "(  7 -8   9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector multiplication assignment";

      initialize();

      RT row2 = row( mat_, 2UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      row2 *= vec;

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  2UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( row2[0] != -4 || row2[1] != 0 || row2[2] != 0 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=  0 || mat_(2,2) != 0 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  0  0  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector multiplication assignment";

      initialize();

      RT row2 = row( mat_, 2UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      row2 *= vec;

      checkSize    ( row2, 4UL );
      checkNonZeros( row2, 1UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( row2[0] != -4 || row2[1] != 0 || row2[2] != 0 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=  0 || mat_(2,2) != 0 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  0  0  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major scalar multiplication assignment";

      initialize();

      RT row2 = row( mat_, 2UL );

      row2 *= 3;

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  2UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( row2[0] != -6 || row2[1] != 0 || row2[2] != -9 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -6 || mat_(2,1) !=  0 || mat_(2,2) != -9 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -6  0 -9  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major SparseRow multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseRow multiplication assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );
      row2 *= row( tmat_, 3UL );

      checkSize    ( row2 , 4UL );
      checkNonZeros( row2 , 1UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( row2[0] != 0 || row2[1] != 0 || row2[2] != -15 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( 0 0 -15 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=   0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -15 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  4 || tmat_(3,2) !=   5 || tmat_(3,3) != -6 ||
          tmat_(4,0) != 7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0   0  0 )\n"
                                     "(  0  1   0  0 )\n"
                                     "(  0  0 -15  0 )\n"
                                     "(  0  4   5 -6 )\n"
                                     "(  7 -8   9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector multiplication assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      row2 *= vec;

      checkSize    ( row2 , 4UL );
      checkNonZeros( row2 , 1UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( row2[0] != -4 || row2[1] != 0 || row2[2] != 0 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) != 0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -4 || tmat_(2,1) !=  0 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) != 5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) != 9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  0  0  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector multiplication assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      row2 *= vec;

      checkSize    ( row2 , 4UL );
      checkNonZeros( row2 , 1UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( row2[0] != -4 || row2[1] != 0 || row2[2] != 0 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) != 0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -4 || tmat_(2,1) !=  0 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) != 5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) != 9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  0  0  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major scalar multiplication assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      row2 *= 3;

      checkSize    ( row2 ,  4UL );
      checkNonZeros( row2 ,  2UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( row2[0] != -6 || row2[1] != 0 || row2[2] != -9 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -6 || tmat_(2,1) !=  0 || tmat_(2,2) != -9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -6  0 -9  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseRow division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the SparseRow class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testDivAssign()
{
   //=====================================================================================
   // Row-major scalar division assignment
   //=====================================================================================

   {
      test_ = "Row-major scalar division assignment";

      initialize();

      RT row2 = row( mat_, 2UL );

      row2 /= 0.5;

      checkSize    ( row2,  4UL );
      checkNonZeros( row2,  2UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( row2[0] != -4 || row2[1] != 0 || row2[2] != -6 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=  0 || mat_(2,2) != -6 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  0 -6  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar division assignment
   //=====================================================================================

   {
      test_ = "Column-major scalar division assignment";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      row2 /= 0.5;

      checkSize    ( row2 ,  4UL );
      checkNonZeros( row2 ,  2UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( row2[0] != -4 || row2[1] != 0 || row2[2] != -6 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -4 || tmat_(2,1) !=  0 || tmat_(2,2) != -6 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  0 -6  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseRow subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the SparseRow class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSubscript()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::operator[]";

      initialize();

      RT row2 = row( mat_, 2UL );

      // Writing the first element
      row2[1] = 9;

      checkSize    ( row2, 4UL );
      checkNonZeros( row2, 3UL );

      if( row2[0] != -2 || row2[1] != 9 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 9 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  9 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the second element
      row2[2] = 0;

      checkSize    ( row2, 4UL );
      checkNonZeros( row2, 2UL );

      if( row2[0] != -2 || row2[1] != 9 || row2[2] != 0 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 9 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != 0 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  9  0  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the third element
      row2[3] = -8;

      checkSize    ( row2, 4UL );
      checkNonZeros( row2, 3UL );

      if( row2[0] != -2 || row2[1] != 9 || row2[2] != 0 || row2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != 0 || mat_(2,3) != -8 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  9  0 -8 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::operator[]";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      // Writing the first element
      row2[1] = 9;

      checkSize    ( row2, 4UL );
      checkNonZeros( row2, 3UL );

      if( row2[0] != -2 || row2[1] != 9 || row2[2] != -3 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 9 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  9 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  9 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the second element
      row2[2] = 0;

      checkSize    ( row2, 4UL );
      checkNonZeros( row2, 2UL );

      if( row2[0] != -2 || row2[1] != 9 || row2[2] != 0 || row2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 9 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) != 0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  9 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) != 5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) != 9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  9  0  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the third element
      row2[3] = -8;

      checkSize    ( row2, 4UL );
      checkNonZeros( row2, 3UL );

      if( row2[0] != -2 || row2[1] != 9 || row2[2] != 0 || row2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row2 << "\n"
             << "   Expected result:\n( -2 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) != 0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  9 || tmat_(2,2) != 0 || tmat_(2,3) != -8 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) != 5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) != 9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  9  0 -8 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseRow iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the SparseRow class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      initialize();

      // Counting the number of elements in 0th column
      {
         test_ = "Row-major iterator subtraction";

         RT row0 = row( mat_, 0UL );
         const size_t number( row0.end() - row0.begin() );

         if( number != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st column
      {
         test_ = "Row-major iterator subtraction";

         RT row1 = row( mat_, 1UL );
         const size_t number( row1.end() - row1.begin() );

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

      // Counting the number of elements in 2nd column
      {
         test_ = "Row-major iterator subtraction";

         RT row2 = row( mat_, 2UL );
         const size_t number( row2.end() - row2.begin() );

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

      // Counting the number of elements in 3rd column
      {
         test_ = "Row-major iterator subtraction";

         RT row3 = row( mat_, 3UL );
         const size_t number( row3.end() - row3.begin() );

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

      // Counting the number of elements in 4th column
      {
         test_ = "Row-major iterator subtraction";

         RT row4 = row( mat_, 4UL );
         const size_t number( row4.end() - row4.begin() );

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
         test_ = "Row-major read-only access via ConstIterator";

         RT row2 = row( mat_, 2UL );
         RT::ConstIterator it ( row2.cbegin() );
         RT::ConstIterator end( row2.cend() );

         if( it == end || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != -3 ) {
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
         test_ = "Row-major assignment via Iterator";

         RT row4 = row( mat_, 4UL );
         int value = 6;

         for( RT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it = value++;
         }

         if( row4[0] != 6 || row4[1] != 7 || row4[2] != 8 || row4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) != 4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  6 || mat_(4,1) != 7 || mat_(4,2) !=  8 || mat_(4,3) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  6  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         RT row4 = row( mat_, 4UL );
         int value = 2;

         for( RT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it += value++;
         }

         if( row4[0] != 8 || row4[1] != 10 || row4[2] != 12 || row4[3] != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  8 || mat_(4,1) != 10 || mat_(4,2) != 12 || mat_(4,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         RT row4 = row( mat_, 4UL );
         int value = 2;

         for( RT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it -= value++;
         }

         if( row4[0] != 6 || row4[1] != 7 || row4[2] != 8 || row4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) != 4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  6 || mat_(4,1) != 7 || mat_(4,2) !=  8 || mat_(4,3) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  6  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         RT row4 = row( mat_, 4UL );
         int value = 1;

         for( RT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it *= value++;
         }

         if( row4[0] != 6 || row4[1] != 14 || row4[2] != 24 || row4[3] != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  6 || mat_(4,1) != 14 || mat_(4,2) != 24 || mat_(4,3) != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         RT row4 = row( mat_, 4UL );

         for( RT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it /= 2;
         }

         if( row4[0] != 3 || row4[1] != 7 || row4[2] != 12 || row4[3] != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 3 7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) != 4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  3 || mat_(4,1) != 7 || mat_(4,2) != 12 || mat_(4,3) != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  3  7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      initialize();

      // Counting the number of elements in 0th column
      {
         test_ = "Row-major iterator subtraction";

         TRT row0 = row( tmat_, 0UL );
         const size_t number( row0.end() - row0.begin() );

         if( number != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st column
      {
         test_ = "Row-major iterator subtraction";

         TRT row1 = row( tmat_, 1UL );
         const size_t number( row1.end() - row1.begin() );

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

      // Counting the number of elements in 2nd column
      {
         test_ = "Row-major iterator subtraction";

         TRT row2 = row( tmat_, 2UL );
         const size_t number( row2.end() - row2.begin() );

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

      // Counting the number of elements in 3rd column
      {
         test_ = "Row-major iterator subtraction";

         TRT row3 = row( tmat_, 3UL );
         const size_t number( row3.end() - row3.begin() );

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

      // Counting the number of elements in 4th column
      {
         test_ = "Row-major iterator subtraction";

         TRT row4 = row( tmat_, 4UL );
         const size_t number( row4.end() - row4.begin() );

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
         test_ = "Column-major read-only access via ConstIterator";

         TRT row2 = row( tmat_, 2UL );
         TRT::ConstIterator it ( row2.cbegin() );
         TRT::ConstIterator end( row2.cend() );

         if( it == end || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != -3 ) {
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
         test_ = "Column-major assignment via Iterator";

         TRT row4 = row( tmat_, 4UL );
         int value = 6;

         for( TRT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it = value++;
         }

         if( row4[0] != 6 || row4[1] != 7 || row4[2] != 8 || row4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != -2 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) != 4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
             tmat_(4,0) !=  6 || tmat_(4,1) != 7 || tmat_(4,2) !=  8 || tmat_(4,3) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  6  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         TRT row4 = row( tmat_, 4UL );
         int value = 2;

         for( TRT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it += value++;
         }

         if( row4[0] != 8 || row4[1] != 10 || row4[2] != 12 || row4[3] != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
             tmat_(4,0) !=  8 || tmat_(4,1) != 10 || tmat_(4,2) != 12 || tmat_(4,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         TRT row4 = row( tmat_, 4UL );
         int value = 2;

         for( TRT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it -= value++;
         }

         if( row4[0] != 6 || row4[1] != 7 || row4[2] != 8 || row4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
             tmat_(4,0) !=  6 || tmat_(4,1) !=  7 || tmat_(4,2) !=  8 || tmat_(4,3) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  6  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         TRT row4 = row( tmat_, 4UL );
         int value = 1;

         for( TRT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it *= value++;
         }

         if( row4[0] != 6 || row4[1] != 14 || row4[2] != 24 || row4[3] != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
             tmat_(4,0) !=  6 || tmat_(4,1) != 14 || tmat_(4,2) != 24 || tmat_(4,3) != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         TRT row4 = row( tmat_, 4UL );

         for( TRT::Iterator it=row4.begin(); it!=row4.end(); ++it ) {
            *it /= 2;
         }

         if( row4[0] != 3 || row4[1] != 7 || row4[2] != 12 || row4[3] != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 3 7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != -2 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) != 4 || tmat_(3,2) !=  5 || tmat_(3,3) != -6 ||
             tmat_(4,0) !=  3 || tmat_(4,1) != 7 || tmat_(4,2) != 12 || tmat_(4,3) != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  3  7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the nonZeros member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of SparseRow. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::nonZeros()";

      initialize();

      // Initialization check
      RT row3 = row( mat_, 3UL );

      checkSize    ( row3, 4UL );
      checkNonZeros( row3, 3UL );

      if( row3[0] != 0 || row3[1] != 4 || row3[2] != 5 || row3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << row3 << "\n"
             << "   Expected result:\n( 0 4 5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse row
      row3[2] = 0;

      checkSize    ( row3, 4UL );
      checkNonZeros( row3, 2UL );

      if( row3[0] != 0 || row3[1] != 4 || row3[2] != 0 || row3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row3 << "\n"
             << "   Expected result:\n( 0 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      mat_(3,0) = 5;

      checkSize    ( row3, 4UL );
      checkNonZeros( row3, 3UL );

      if( row3[0] != 5 || row3[1] != 4 || row3[2] != 0 || row3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row3 << "\n"
             << "   Expected result:\n( 5 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::nonZeros()";

      initialize();

      // Initialization check
      TRT row3 = row( tmat_, 3UL );

      checkSize    ( row3, 4UL );
      checkNonZeros( row3, 3UL );

      if( row3[0] != 0 || row3[1] != 4 || row3[2] != 5 || row3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << row3 << "\n"
             << "   Expected result:\n( 0 4 5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse row
      row3[2] = 0;

      checkSize    ( row3, 4UL );
      checkNonZeros( row3, 2UL );

      if( row3[0] != 0 || row3[1] != 4 || row3[2] != 0 || row3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row3 << "\n"
             << "   Expected result:\n( 0 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      tmat_(3,0) = 5;

      checkSize    ( row3, 4UL );
      checkNonZeros( row3, 3UL );

      if( row3[0] != 5 || row3[1] != 4 || row3[2] != 0 || row3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << row3 << "\n"
             << "   Expected result:\n( 5 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reset member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of SparseRow. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::reset()";

      initialize();

      // Resetting the 0th row
      {
         RT row0 = row( mat_, 0UL );
         row0.reset();

         checkSize    ( row0,  4UL );
         checkNonZeros( row0,  0UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 10UL );

         if( row0[0] != 0 || row0[1] != 0 || row0[2] != 0 || row0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 0th row failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st row
      {
         RT row1 = row( mat_, 1UL );
         row1.reset();

         checkSize    ( row1, 4UL );
         checkNonZeros( row1, 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 9UL );

         if( row1[0] != 0 || row1[1] != 0 || row1[2] != 0 || row1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st row failed\n"
                << " Details:\n"
                << "   Result:\n" << row1 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd row
      {
         RT row2 = row( mat_, 2UL );
         row2.reset();

         checkSize    ( row2, 4UL );
         checkNonZeros( row2, 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 7UL );

         if( row2[0] != 0 || row2[1] != 0 || row2[2] != 0 || row2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd row failed\n"
                << " Details:\n"
                << "   Result:\n" << row2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 3rd row
      {
         RT row3 = row( mat_, 3UL );
         row3.reset();

         checkSize    ( row3, 4UL );
         checkNonZeros( row3, 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 4UL );

         if( row3[0] != 0 || row3[1] != 0 || row3[2] != 0 || row3[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 4th row
      {
         RT row4 = row( mat_, 4UL );
         row4.reset();

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 0UL );

         if( row4[0] != 0 || row4[1] != 0 || row4[2] != 0 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 4th row failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::reset()";

      initialize();

      // Resetting the 0th row
      {
         TRT row0 = row( tmat_, 0UL );
         row0.reset();

         checkSize    ( row0 ,  4UL );
         checkNonZeros( row0 ,  0UL );
         checkRows    ( tmat_,  5UL );
         checkColumns ( tmat_,  4UL );
         checkNonZeros( tmat_, 10UL );

         if( row0[0] != 0 || row0[1] != 0 || row0[2] != 0 || row0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 0th row failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st row
      {
         TRT row1 = row( tmat_, 1UL );
         row1.reset();

         checkSize    ( row1 , 4UL );
         checkNonZeros( row1 , 0UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 9UL );

         if( row1[0] != 0 || row1[1] != 0 || row1[2] != 0 || row1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st row failed\n"
                << " Details:\n"
                << "   Result:\n" << row1 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd row
      {
         TRT row2 = row( tmat_, 2UL );
         row2.reset();

         checkSize    ( row2 , 4UL );
         checkNonZeros( row2 , 0UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 7UL );

         if( row2[0] != 0 || row2[1] != 0 || row2[2] != 0 || row2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd row failed\n"
                << " Details:\n"
                << "   Result:\n" << row2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 3rd row
      {
         TRT row3 = row( tmat_, 3UL );
         row3.reset();

         checkSize    ( row3 , 4UL );
         checkNonZeros( row3 , 0UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 4UL );

         if( row3[0] != 0 || row3[1] != 0 || row3[2] != 0 || row3[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 4th row
      {
         TRT row4 = row( tmat_, 4UL );
         row4.reset();

         checkSize    ( row4 , 4UL );
         checkNonZeros( row4 , 0UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 0UL );

         if( row4[0] != 0 || row4[1] != 0 || row4[2] != 0 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 4th row failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the append member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the append member function of SparseRow. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::append()";

      MT mat( 3UL, 9UL );

      RT row1 = row( mat, 1UL );
      row1.reserve( 4UL );

      // Appending one non-zero element
      row1.append( 1UL, 1 );

      checkSize    ( row1, 9UL );
      checkCapacity( row1, 4UL );
      checkNonZeros( row1, 1UL );

      if( row1[1] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 1 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Appending three more non-zero elements
      row1.append( 3UL, 2 );
      row1.append( 4UL, 3 );
      row1.append( 8UL, 4 );

      checkSize    ( row1, 9UL );
      checkCapacity( row1, 4UL );
      checkNonZeros( row1, 4UL );

      if( row1[1] != 1 || row1[3] != 2 || row1[4] != 3 || row1[8] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 1 0 2 3 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::append()";

      TMT mat( 3UL, 9UL );

      TRT row1 = row( mat, 1UL );
      row1.reserve( 4UL );

      // Appending one non-zero element
      row1.append( 1UL, 1 );

      checkSize    ( row1, 9UL );
      checkCapacity( row1, 4UL );
      checkNonZeros( row1, 1UL );

      if( row1[1] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 1 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Appending three more non-zero elements
      row1.append( 3UL, 2 );
      row1.append( 4UL, 3 );
      row1.append( 8UL, 4 );

      checkSize    ( row1, 9UL );
      checkCapacity( row1, 4UL );
      checkNonZeros( row1, 4UL );

      if( row1[1] != 1 || row1[3] != 2 || row1[4] != 3 || row1[8] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 1 0 2 3 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the insert member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the insert member function of SparseRow. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::insert()";

      initialize();

      RT row0 = row( mat_, 0UL );

      // Inserting a non-zero element at the end of the row
      {
         RT::Iterator pos = row0.insert( 3UL, 1 );

         checkSize    ( row0,  4UL );
         checkNonZeros( row0,  1UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 11UL );

         if( pos->value() != 1 || pos->index() != 3UL ) {
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

         if( row0[0] != 0 || row0[1] != 0 || row0[2] != 0 || row0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the beginning of the row
      {
         RT::Iterator pos = row0.insert( 0UL, 2 );

         checkSize    ( row0,  4UL );
         checkNonZeros( row0,  2UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 12UL );

         if( pos->value() != 2 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( row0[0] != 2 || row0[1] != 0 || row0[2] != 0 || row0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 2 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the center of the row
      {
         RT::Iterator pos = row0.insert( 2UL, 3 );

         checkSize    ( row0,  4UL );
         checkNonZeros( row0,  3UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 13UL );

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

         if( row0[0] != 2 || row0[1] != 0 || row0[2] != 3 || row0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 2 0 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         row0.insert( 3UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << row0 << "\n"
             << "   Expected result:\n( 2 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::insert()";

      initialize();

      TRT row0 = row( tmat_, 0UL );

      // Inserting a non-zero element at the end of the row
      {
         TRT::Iterator pos = row0.insert( 3UL, 1 );

         checkSize    ( row0 ,  4UL );
         checkNonZeros( row0 ,  1UL );
         checkRows    ( tmat_,  5UL );
         checkColumns ( tmat_,  4UL );
         checkNonZeros( tmat_, 11UL );

         if( pos->value() != 1 || pos->index() != 3UL ) {
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

         if( row0[0] != 0 || row0[1] != 0 || row0[2] != 0 || row0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the beginning of the row
      {
         TRT::Iterator pos = row0.insert( 0UL, 2 );

         checkSize    ( row0 ,  4UL );
         checkNonZeros( row0 ,  2UL );
         checkRows    ( tmat_,  5UL );
         checkColumns ( tmat_,  4UL );
         checkNonZeros( tmat_, 12UL );

         if( pos->value() != 2 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( row0[0] != 2 || row0[1] != 0 || row0[2] != 0 || row0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 2 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the center of the row
      {
         TRT::Iterator pos = row0.insert( 2UL, 3 );

         checkSize    ( row0 ,  4UL );
         checkNonZeros( row0 ,  3UL );
         checkRows    ( tmat_,  5UL );
         checkColumns ( tmat_,  4UL );
         checkNonZeros( tmat_, 13UL );

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

         if( row0[0] != 2 || row0[1] != 0 || row0[2] != 3 || row0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row0 << "\n"
                << "   Expected result:\n( 2 0 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         row0.insert( 3UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << row0 << "\n"
             << "   Expected result:\n( 2 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the erase member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the erase member function of SparseRow. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::erase( size_t )";

      initialize();

      RT row4 = row( mat_, 4UL );

      // Erasing the non-zero element at the end of the row
      row4.erase( 3UL );

      checkSize    ( row4, 4UL );
      checkNonZeros( row4, 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( row4[0] != 7 || row4[1] != -8 || row4[2] != 9 || row4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 7 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the row
      row4.erase( size_t(0) );

      checkSize    ( row4, 4UL );
      checkNonZeros( row4, 2UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 8UL );

      if( row4[0] != 0 || row4[1] != -8 || row4[2] != 9 || row4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the row
      row4.erase( 2UL );

      checkSize    ( row4, 4UL );
      checkNonZeros( row4, 1UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 7UL );

      if( row4[0] != 0 || row4[1] != -8 || row4[2] != 0 || row4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      row4.erase( 3UL );

      checkSize    ( row4, 4UL );
      checkNonZeros( row4, 1UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 7UL );

      if( row4[0] != 0 || row4[1] != -8 || row4[2] != 0 || row4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::erase( Iterator )";

      initialize();

      RT row4 = row( mat_, 4UL );

      // Erasing the non-zero element at the end of the row
      {
         RT::Iterator pos = row4.erase( row4.find( 3UL ) );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 3UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 9UL );

         if( pos != row4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 7 || row4[1] != -8 || row4[2] != 9 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 7 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the row
      {
         RT::Iterator pos = row4.erase( row4.find( 0UL ) );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 2UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 8UL );

         if( pos->value() != -8 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -8\n"
                << "   Expected index:  1\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 0 || row4[1] != -8 || row4[2] != 9 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the row
      {
         RT::Iterator pos = row4.erase( row4.find( 2UL ) );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 1UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 7UL );

         if( pos != row4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 0 || row4[1] != -8 || row4[2] != 0 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         RT::Iterator pos = row4.erase( row4.find( 3UL ) );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 1UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 7UL );

         if( pos != row4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 0 || row4[1] != -8 || row4[2] != 0 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::erase( Iterator, Iterator )";

      initialize();

      // Erasing the 2nd row
      {
         RT row2 = row( mat_, 2UL );

         RT::Iterator pos = row2.erase( row2.begin(), row2.end() );

         checkSize    ( row2, 4UL );
         checkNonZeros( row2, 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 8UL );

         if( pos != row2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row2[0] != 0 || row2[1] != -0 || row2[2] != 0 || row2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the row failed\n"
                << " Details:\n"
                << "   Result:\n" << row2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 4th row
      {
         RT row4 = row( mat_, 4UL );

         RT::Iterator pos = row4.erase( row4.begin(), row4.find( 2UL ) );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 2UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 6UL );

         if( pos->value() != 9 || pos->index() != 2 ) {
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

         if( row4[0] != 0 || row4[1] != 0 || row4[2] != 9 || row4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial row failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 0 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 4th row
      {
         RT row4 = row( mat_, 4UL );

         RT::Iterator pos = row4.erase( row4.find( 2UL ), row4.end() );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 4UL );

         if( pos != row4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 0 || row4[1] != 0 || row4[2] != 0 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial row failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         RT row3 = row( mat_, 3UL );

         RT::Iterator pos = row3.erase( row3.find( 1UL ), row3.find( 1UL ) );

         checkSize    ( row3, 4UL );
         checkNonZeros( row3, 3UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 4UL );

         if( pos != row3.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row3[0] != 0 || row3[1] != 4 || row3[2] != 5 || row3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::erase( size_t )";

      initialize();

      TRT row4 = row( tmat_, 4UL );

      // Erasing the non-zero element at the end of the row
      row4.erase( 3UL );

      checkSize    ( row4 , 4UL );
      checkNonZeros( row4 , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( row4[0] != 7 || row4[1] != -8 || row4[2] != 9 || row4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 7 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the row
      row4.erase( 0UL );

      checkSize    ( row4 , 4UL );
      checkNonZeros( row4 , 2UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 8UL );

      if( row4[0] != 0 || row4[1] != -8 || row4[2] != 9 || row4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the row
      row4.erase( 2UL );

      checkSize    ( row4 , 4UL );
      checkNonZeros( row4 , 1UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 7UL );

      if( row4[0] != 0 || row4[1] != -8 || row4[2] != 0 || row4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      row4.erase( 3UL );

      checkSize    ( row4 , 4UL );
      checkNonZeros( row4 , 1UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 7UL );

      if( row4[0] != 0 || row4[1] != -8 || row4[2] != 0 || row4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << row4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::erase( Iterator )";

      initialize();

      TRT row4 = row( tmat_, 4UL );

      // Erasing the non-zero element at the end of the row
      {
         TRT::Iterator pos = row4.erase( row4.find( 3UL ) );

         checkSize    ( row4 , 4UL );
         checkNonZeros( row4 , 3UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 9UL );

         if( pos != row4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 7 || row4[1] != -8 || row4[2] != 9 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 7 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the row
      {
         TRT::Iterator pos = row4.erase( row4.find( 0UL ) );

         checkSize    ( row4 , 4UL );
         checkNonZeros( row4 , 2UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 8UL );

         if( pos->value() != -8 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -8\n"
                << "   Expected index:  1\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 0 || row4[1] != -8 || row4[2] != 9 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the row
      {
         TRT::Iterator pos = row4.erase( row4.find( 2UL ) );

         checkSize    ( row4 , 4UL );
         checkNonZeros( row4 , 1UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != row4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 0 || row4[1] != -8 || row4[2] != 0 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         TRT::Iterator pos = row4.erase( row4.find( 3UL ) );

         checkSize    ( row4 , 4UL );
         checkNonZeros( row4 , 1UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != row4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 0 || row4[1] != -8 || row4[2] != 0 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::erase( Iterator, Iterator )";

      initialize();

      // Erasing the 2nd row
      {
         TRT row2 = row( tmat_, 2UL );

         TRT::Iterator pos = row2.erase( row2.begin(), row2.end() );

         checkSize    ( row2 , 4UL );
         checkNonZeros( row2 , 0UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 8UL );

         if( pos != row2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row2[0] != 0 || row2[1] != -0 || row2[2] != 0 || row2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the row failed\n"
                << " Details:\n"
                << "   Result:\n" << row2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 4th row
      {
         TRT row4 = row( tmat_, 4UL );

         TRT::Iterator pos = row4.erase( row4.begin(), row4.find( 2UL ) );

         checkSize    ( row4 , 4UL );
         checkNonZeros( row4 , 2UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 6UL );

         if( pos->value() != 9 || pos->index() != 2 ) {
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

         if( row4[0] != 0 || row4[1] != 0 || row4[2] != 9 || row4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial row failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 0 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 4th row
      {
         TRT row4 = row( tmat_, 4UL );

         TRT::Iterator pos = row4.erase( row4.find( 2UL ), row4.end() );

         checkSize    ( row4 , 4UL );
         checkNonZeros( row4 , 0UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 4UL );

         if( pos != row4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row4[0] != 0 || row4[1] != 0 || row4[2] != 0 || row4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial row failed\n"
                << " Details:\n"
                << "   Result:\n" << row4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         TRT row3 = row( tmat_, 3UL );

         TRT::Iterator pos = row3.erase( row3.find( 1UL ), row3.find( 1UL ) );

         checkSize    ( row3 , 4UL );
         checkNonZeros( row3 , 3UL );
         checkRows    ( tmat_, 5UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 4UL );

         if( pos != row3.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( row3[0] != 0 || row3[1] != 4 || row3[2] != 5 || row3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reserve member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reserve member function of SparseRow. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::reserve()";

      MT mat( 3UL, 20UL );

      RT row0 = row( mat, 0UL );

      // Increasing the capacity of the row
      row0.reserve( 10UL );

      checkSize    ( row0, 20UL );
      checkCapacity( row0, 10UL );
      checkNonZeros( row0,  0UL );

      // Further increasing the capacity of the row
      row0.reserve( 15UL );

      checkSize    ( row0, 20UL );
      checkCapacity( row0, 15UL );
      checkNonZeros( row0,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::reserve()";

      TMT mat( 3UL, 20UL );

      TRT row0 = row( mat, 0UL );

      // Increasing the capacity of the row
      row0.reserve( 10UL );

      checkSize    ( row0, 20UL );
      checkCapacity( row0, 10UL );
      checkNonZeros( row0,  0UL );

      // Further increasing the capacity of the row
      row0.reserve( 15UL );

      checkSize    ( row0, 20UL );
      checkCapacity( row0, 15UL );
      checkNonZeros( row0,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of SparseRow. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScale()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::scale()";

      initialize();

      // Integral scaling the 3rd row
      {
         RT row3 = row( mat_, 3UL );
         row3.scale( 3 );

         checkSize    ( row3,  4UL );
         checkNonZeros( row3,  3UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 10UL );

         if( row3[0] != 0 || row3[1] != 12 || row3[2] != 15 || row3[3] != -18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 12 15 -18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=   0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=   0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=   0 ||
             mat_(3,0) !=  0 || mat_(3,1) != 12 || mat_(3,2) != 15 || mat_(3,3) != -18 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) !=  10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0   0   0   0 )\n"
                                        "(  0   1   0   0 )\n"
                                        "( -2   0  -3   0 )\n"
                                        "(  0  12  15 -18 )\n"
                                        "(  7  -8   9  10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Floating point scaling the 3rd row
      {
         RT row3 = row( mat_, 3UL );
         row3.scale( 0.5 );

         checkSize    ( row3,  4UL );
         checkNonZeros( row3,  3UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 10UL );

         if( row3[0] != 0 || row3[1] != 6 || row3[2] != 7 || row3[3] != -9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 6 7 -9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  6 || mat_(3,2) !=  7 || mat_(3,3) != -9 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0   0   0   0 )\n"
                                        "(  0   1   0   0 )\n"
                                        "( -2   0  -3   0 )\n"
                                        "(  0   6   7  -9 )\n"
                                        "(  7  -8   9  10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::scale()";

      initialize();

      // Integral scaling the 3rd row
      {
         TRT row3 = row( tmat_, 3UL );
         row3.scale( 3 );

         checkSize    ( row3 ,  4UL );
         checkNonZeros( row3 ,  3UL );
         checkRows    ( tmat_,  5UL );
         checkColumns ( tmat_,  4UL );
         checkNonZeros( tmat_, 10UL );

         if( row3[0] != 0 || row3[1] != 12 || row3[2] != 15 || row3[3] != -18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 12 15 -18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=   0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=   0 ||
             tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=   0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) != 12 || tmat_(3,2) != 15 || tmat_(3,3) != -18 ||
             tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) !=  10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0   0   0   0 )\n"
                                        "(  0   1   0   0 )\n"
                                        "( -2   0  -3   0 )\n"
                                        "(  0  12  15 -18 )\n"
                                        "(  7  -8   9  10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Floating point scaling the 3rd row
      {
         TRT row3 = row( tmat_, 3UL );
         row3.scale( 0.5 );

         checkSize    ( row3 ,  4UL );
         checkNonZeros( row3 ,  3UL );
         checkRows    ( tmat_,  5UL );
         checkColumns ( tmat_,  4UL );
         checkNonZeros( tmat_, 10UL );

         if( row3[0] != 0 || row3[1] != 6 || row3[2] != 7 || row3[3] != -9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << row3 << "\n"
                << "   Expected result:\n( 0 6 7 -9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) !=  6 || tmat_(3,2) !=  7 || tmat_(3,3) != -9 ||
             tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0   0   0   0 )\n"
                                        "(  0   1   0   0 )\n"
                                        "( -2   0  -3   0 )\n"
                                        "(  0   6   7  -9 )\n"
                                        "(  7  -8   9  10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the find member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the find member function of SparseRow. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::find()";

      initialize();

      RT row2 = row( mat_, 2UL );

      // Searching for the first element
      {
         RT::Iterator pos = row2.find( 0UL );

         if( pos == row2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 0 || pos->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         RT::Iterator pos = row2.find( 2UL );

         if( pos == row2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         RT::Iterator pos = row2.find( 1UL );

         if( pos != row2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::find()";

      initialize();

      TRT row2 = row( tmat_, 2UL );

      // Searching for the first element
      {
         TRT::Iterator pos = row2.find( 0UL );

         if( pos == row2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 0 || pos->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         TRT::Iterator pos = row2.find( 2UL );

         if( pos == row2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         TRT::Iterator pos = row2.find( 1UL );

         if( pos != row2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the lowerBound member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the lowerBound member function of SparseRow. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::lowerBound()";

      initialize();

      RT row1 = row( mat_, 1UL );

      // Determining the lower bound for index 0
      {
         RT::Iterator pos = row1.lowerBound( 0UL );

         if( pos == row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current row:\n" << row1 << "\n";
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
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 1
      {
         RT::Iterator pos = row1.lowerBound( 1UL );

         if( pos == row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current row:\n" << row1 << "\n";
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
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 2
      {
         RT::Iterator pos = row1.lowerBound( 2UL );

         if( pos != row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::lowerBound()";

      initialize();

      TRT row1 = row( tmat_, 1UL );

      // Determining the lower bound for index 0
      {
         TRT::Iterator pos = row1.lowerBound( 0UL );

         if( pos == row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current row:\n" << row1 << "\n";
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
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 1
      {
         TRT::Iterator pos = row1.lowerBound( 1UL );

         if( pos == row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current row:\n" << row1 << "\n";
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
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 2
      {
         TRT::Iterator pos = row1.lowerBound( 2UL );

         if( pos != row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the upperBound member function of SparseRow.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the upperBound member function of SparseRow. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseRow::upperBound()";

      initialize();

      RT row1 = row( mat_, 1UL );

      // Determining the upper bound for index 0
      {
         RT::Iterator pos = row1.upperBound( 0UL );

         if( pos == row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current row:\n" << row1 << "\n";
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
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 1
      {
         RT::Iterator pos = row1.upperBound( 1UL );

         if( pos != row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 2
      {
         RT::Iterator pos = row1.upperBound( 2UL );

         if( pos != row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseRow::upperBound()";

      initialize();

      TRT row1 = row( tmat_, 1UL );

      // Determining the upper bound for index 0
      {
         TRT::Iterator pos = row1.upperBound( 0UL );

         if( pos == row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current row:\n" << row1 << "\n";
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
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 1
      {
         TRT::Iterator pos = row1.upperBound( 1UL );

         if( pos != row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 2
      {
         TRT::Iterator pos = row1.upperBound( 2UL );

         if( pos != row1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDefault function with the SparseRow class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDefault function with the SparseRow class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      initialize();

      // isDefault with default row
      {
         RT row0 = row( mat_, 0UL );

         if( isDefault( row0 ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default row
      {
         RT row1 = row( mat_, 1UL );

         if( isDefault( row1 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDefault() function";

      initialize();

      // isDefault with default row
      {
         TRT row0 = row( tmat_, 0UL );

         if( isDefault( row0 ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default row
      {
         TRT row1 = row( tmat_, 1UL );

         if( isDefault( row1 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isnan function with the SparseRow class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isnan function with the SparseRow class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsNan()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isnan() function";

      typedef blaze::CompressedMatrix<float,blaze::rowMajor>  MatrixType;
      typedef blaze::SparseRow<MatrixType>                    RowType;

      MatrixType mat( mat_ );

      checkRows    ( mat,  5UL );
      checkColumns ( mat,  4UL );
      checkNonZeros( mat, 10UL );

      // isnan with empty row
      {
         RowType row0 = row( mat, 0UL );

         checkSize    ( row0, 4UL );
         checkNonZeros( row0, 0UL );

         if( blaze::isnan( row0 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with partially filled row
      {
         RowType row2 = row( mat, 2UL );

         checkSize    ( row2, 4UL );
         checkNonZeros( row2, 2UL );

         if( blaze::isnan( row2 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with fully filled row
      {
         RowType row4 = row( mat, 4UL );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 4UL );

         if( blaze::isnan( row4 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row4 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isnan() function";

      typedef blaze::CompressedMatrix<float,blaze::columnMajor>  MatrixType;
      typedef blaze::SparseRow<MatrixType>                       RowType;

      MatrixType mat( mat_ );

      checkRows    ( mat,  5UL );
      checkColumns ( mat,  4UL );
      checkNonZeros( mat, 10UL );

      // isnan with empty row
      {
         RowType row0 = row( mat, 0UL );

         checkSize    ( row0, 4UL );
         checkNonZeros( row0, 0UL );

         if( blaze::isnan( row0 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with partially filled row
      {
         RowType row2 = row( mat, 2UL );

         checkSize    ( row2, 4UL );
         checkNonZeros( row2, 2UL );

         if( blaze::isnan( row2 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with fully filled row
      {
         RowType row4 = row( mat, 4UL );

         checkSize    ( row4, 4UL );
         checkNonZeros( row4, 4UL );

         if( blaze::isnan( row4 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Row:\n" << row4 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the min function with the SparseRow class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the min function used with the SparseRow class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMinimum()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major min() function";

      initialize();

      // Computing the minimum of the 0th row
      {
         const int minimum = min( row( mat_, 0UL ) );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 0th row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 1st row
      {
         const int minimum = min( row( mat_, 1UL ) );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 1st row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 2nd row
      {
         const int minimum = min( row( mat_, 2UL ) );

         if( minimum != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 2nd row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 3rd row
      {
         const int minimum = min( row( mat_, 3UL ) );

         if( minimum != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 3rd row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 4th row
      {
         const int minimum = min( row( mat_, 4UL ) );

         if( minimum != -8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 4th row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -8\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major min() function";

      initialize();

      // Computing the minimum of the 0th row
      {
         const int minimum = min( row( tmat_, 0UL ) );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 0th row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 1st row
      {
         const int minimum = min( row( tmat_, 1UL ) );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 1st row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 2nd row
      {
         const int minimum = min( row( tmat_, 2UL ) );

         if( minimum != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 2nd row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 3rd row
      {
         const int minimum = min( row( tmat_, 3UL ) );

         if( minimum != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 3rd row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 4th row
      {
         const int minimum = min( row( tmat_, 4UL ) );

         if( minimum != -8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 4th row failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -8\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the max function with the SparseRow class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the max function used with the SparseRow class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMaximum()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major max() function";

      initialize();

      // Computing the maximum of the 0th row
      {
         const int maximum = max( row( mat_, 0UL ) );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 0th row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 1st row
      {
         const int maximum = max( row( mat_, 1UL ) );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 1st row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 2nd row
      {
         const int maximum = max( row( mat_, 2UL ) );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 2nd row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 3rd row
      {
         const int maximum = max( row( mat_, 3UL ) );

         if( maximum != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 3rd row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 4th row
      {
         const int maximum = max( row( mat_, 4UL ) );

         if( maximum != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 4th row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 10\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major max() function";

      initialize();

      // Computing the maximum of the 0th row
      {
         const int maximum = max( row( tmat_, 0UL ) );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 0th row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 1st row
      {
         const int maximum = max( row( tmat_, 1UL ) );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 1st row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 2nd row
      {
         const int maximum = max( row( tmat_, 2UL ) );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 2nd row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 3rd row
      {
         const int maximum = max( row( tmat_, 3UL ) );

         if( maximum != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 3rd row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 4th row
      {
         const int maximum = max( row( tmat_, 4UL ) );

         if( maximum != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 4th row failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 10\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subvector function with the SparseRow class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subvector function used with the SparseRow class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubvector()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major subvector() function";

      initialize();

      typedef blaze::SparseSubvector<RT>  SubvectorType;

      RT row1 = row( mat_, 1UL );
      SubvectorType sv = subvector( row1, 0UL, 4UL );

      if( sv[1] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << sv[1] << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv.begin()->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << sv.begin()->value() << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major subvector() function";

      initialize();

      typedef blaze::SparseSubvector<TRT>  SubvectorType;

      TRT row1 = row( tmat_, 1UL );
      SubvectorType sv = subvector( row1, 0UL, 4UL );

      if( sv[1] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << sv[1] << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv.begin()->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << sv.begin()->value() << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member matrices to specific predetermined values.
*/
void ClassTest::initialize()
{
   // Initializing the row-major compressed matrix
   mat_.reset();
   mat_(1,1) =  1;
   mat_(2,0) = -2;
   mat_(2,2) = -3;
   mat_(3,1) =  4;
   mat_(3,2) =  5;
   mat_(3,3) = -6;
   mat_(4,0) =  7;
   mat_(4,1) = -8;
   mat_(4,2) =  9;
   mat_(4,3) = 10;

   // Initializing the column-major compressed matrix
   tmat_.reset();
   tmat_(1,1) =  1;
   tmat_(2,0) = -2;
   tmat_(2,2) = -3;
   tmat_(3,1) =  4;
   tmat_(3,2) =  5;
   tmat_(3,3) = -6;
   tmat_(4,0) =  7;
   tmat_(4,1) = -8;
   tmat_(4,2) =  9;
   tmat_(4,3) = 10;
}
//*************************************************************************************************

} // namespace sparserow

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
   std::cout << "   Running SparseRow class test..." << std::endl;

   try
   {
      RUN_SPARSEROW_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SparseRow class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
