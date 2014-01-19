//=================================================================================================
/*!
//  \file src/mathtest/sparsecolumn/ClassTest.cpp
//  \brief Source file for the SparseColumn class test
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
#include <blazetest/mathtest/sparsecolumn/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace sparsecolumn {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SparseColumn class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
   : mat_ ( 4UL, 5UL )
   , tmat_( 4UL, 5UL )
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
/*!\brief Test of the SparseColumn constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the SparseColumn class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn constructor";

      initialize();

      // 0th matrix column
      {
         CT col0 = column( mat_, 0UL );

         checkSize    ( col0, 4UL );
         checkNonZeros( col0, 0UL );

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 0th sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st matrix column
      {
         CT col1 = column( mat_, 1UL );

         checkSize    ( col1, 4UL );
         checkNonZeros( col1, 1UL );

         if( col1[0] != 0 || col1[1] != 1 || col1[2] != 0 || col1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col1 << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd matrix column
      {
         CT col2 = column( mat_, 2UL );

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 2UL );

         if( col2[0] != -2 || col2[1] != 0 || col2[2] != -3 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( -2 0 -3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd matrix column
      {
         CT col3 = column( mat_, 3UL );

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 3UL );

         if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th matrix column
      {
         CT col4 = column( mat_, 4UL );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 4UL );

         if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 4th sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 7 -8 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn constructor";

      initialize();

      // 0th matrix column
      {
         TCT col0 = column( tmat_, 0UL );

         checkSize    ( col0, 4UL );
         checkNonZeros( col0, 0UL );

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 0th sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st matrix column
      {
         TCT col1 = column( tmat_, 1UL );

         checkSize    ( col1, 4UL );
         checkNonZeros( col1, 1UL );

         if( col1[0] != 0 || col1[1] != 1 || col1[2] != 0 || col1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col1 << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd matrix column
      {
         TCT col2 = column( tmat_, 2UL );

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 2UL );

         if( col2[0] != -2 || col2[1] != 0 || col2[2] != -3 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( -2 0 -3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd matrix column
      {
         TCT col3 = column( tmat_, 3UL );

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 3UL );

         if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th matrix column
      {
         TCT col4 = column( tmat_, 4UL );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 4UL );

         if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 4th sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 7 -8 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseColumn assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the SparseColumn class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn copy assignment";

      initialize();

      CT col1 = column( mat_, 1UL );
      col1 = column( mat_, 2UL );

      checkSize    ( col1,  4UL );
      checkNonZeros( col1,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col1[0] != -2 || col1[1] != 0 || col1[2] != -3 || col1[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( -2 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != -2 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) !=  0 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != -3 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0 -2 -2  0  7 )\n"
                                     "( 0  0  0  4 -8 )\n"
                                     "( 0 -3 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector assignment";

      initialize();

      CT col1 = column( mat_, 1UL );

      blaze::DynamicVector<int,blaze::columnVector> vec1( 4UL, 0 );
      vec1[1] = 8;
      vec1[3] = 9;

      col1 = vec1;

      checkSize    ( col1,  4UL );
      checkNonZeros( col1,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 || col1[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 8 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 9 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  8  0  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  9  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector assignment";

      initialize();

      CT col4 = column( mat_, 4UL );

      blaze::CompressedVector<int,blaze::columnVector> vec1( 4UL );
      vec1[3] = 9;

      col4 = vec1;

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 7UL );

      if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) != 0 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 0 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 0 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  0 )\n"
                                     "( 0  1  0  4  0 )\n"
                                     "( 0  0 -3  5  0 )\n"
                                     "( 0  0  0 -6  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn copy assignment";

      initialize();

      TCT col1 = column( tmat_, 1UL );
      col1 = column( tmat_, 2UL );

      checkSize    ( col1 ,  4UL );
      checkNonZeros( col1 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col1[0] != -2 || col1[1] != 0 || col1[2] != -3 || col1[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( -2 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != -2 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != -3 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 -2 -2  0  7 )\n"
                                     "( 0  0  0  4 -8 )\n"
                                     "( 0 -3 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector assignment";

      initialize();

      TCT col1 = column( tmat_, 1UL );

      blaze::DynamicVector<int,blaze::columnVector> vec1( 4UL, 0 );
      vec1[1] = 8;
      vec1[3] = 9;

      col1 = vec1;

      checkSize    ( col1 ,  4UL );
      checkNonZeros( col1 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 || col1[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 8 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 9 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  8  0  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  9  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector assignment";

      initialize();

      TCT col4 = column( tmat_, 4UL );

      blaze::CompressedVector<int,blaze::columnVector> vec1( 4UL );
      vec1[3] = 9;

      col4 = vec1;

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) != 0 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 0 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  0 )\n"
                                     "( 0  1  0  4  0 )\n"
                                     "( 0  0 -3  5  0 )\n"
                                     "( 0  0  0 -6  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseColumn addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the SparseColumn class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Row-major SparseColumn addition assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn addition assignment";

      initialize();

      CT col2 = column( mat_, 2UL );
      col2 += column( mat_, 3UL );

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( col2[0] != -2 || col2[1] != 4 || col2[2] != 2 || col2[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 4 2 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  2 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != -6 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0  2  5  9 )\n"
                                     "( 0  0 -6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector addition assignment";

      initialize();

      CT col2 = column( mat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      col2 += vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != 0 || col2[1] != -4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) != -4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector addition assignment";

      initialize();

      CT col2 = column( mat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 += vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != 0 || col2[1] != -4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) != -4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major SparseColumn addition assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn addition assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );
      col2 += column( tmat_, 3UL );

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( col2[0] != -2 || col2[1] != 4 || col2[2] != 2 || col2[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 4 2 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  2 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != -6 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0  2  5  9 )\n"
                                     "( 0  0 -6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector addition assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      col2 += vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != 0 || col2[1] != -4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) != -4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector addition assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 += vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  3UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col2[0] != 0 || col2[1] != -4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) != -4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseColumn subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the SparseColumn class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Row-major SparseColumn subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn subtraction assignment";

      initialize();

      CT col2 = column( mat_, 2UL );
      col2 -= column( mat_, 3UL );

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( col2[0] != -2 || col2[1] != -4 || col2[2] != -8 || col2[3] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 -4 -8 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) != -4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -8 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  6 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -8  5  9 )\n"
                                     "( 0  0  6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector subtraction assignment";

      initialize();

      CT col2 = column( mat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      col2 -= vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  3UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col2[0] != -4 || col2[1] != 4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector subtraction assignment";

      initialize();

      CT col2 = column( mat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 -= vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  3UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col2[0] != -4 || col2[1] != 4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major SparseColumn subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn subtraction assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );
      col2 -= column( tmat_, 3UL );

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( col2[0] != -2 || col2[1] != -4 || col2[2] != -8 || col2[3] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 -4 -8 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) != -4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -8 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  6 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -8  5  9 )\n"
                                     "( 0  0  6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector subtraction assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      col2 -= vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  3UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col2[0] != -4 || col2[1] != 4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector subtraction assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 -= vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  3UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col2[0] != -4 || col2[1] != 4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseColumn multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the SparseColumn class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Row-major SparseColumn multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn multiplication assignment";

      initialize();

      CT col2 = column( mat_, 2UL );
      col2 *= column( mat_, 3UL );

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col2[0] != 0 || col2[1] != 0 || col2[2] != -15 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 0 -15 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=   0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -15 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=   0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0   0  0  7 )\n"
                                     "( 0  1   0  4 -8 )\n"
                                     "( 0  0 -15  5  9 )\n"
                                     "( 0  0   0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector multiplication assignment";

      initialize();

      CT col2 = column( mat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      col2 *= vec;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector multiplication assignment";

      initialize();

      CT col2 = column( mat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 *= vec;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major scalar multiplication assignment";

      initialize();

      CT col2 = column( mat_, 2UL );

      col2 *= 3;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != -6 || col2[1] != 0 || col2[2] != -9 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -6 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -9 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major SparseColumn multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn multiplication assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );
      col2 *= column( tmat_, 3UL );

      checkSize    ( col2 , 4UL );
      checkNonZeros( col2 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( col2[0] != 0 || col2[1] != 0 || col2[2] != -15 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 0 -15 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=   0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -15 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=   0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0   0  0  7 )\n"
                                     "( 0  1   0  4 -8 )\n"
                                     "( 0  0 -15  5  9 )\n"
                                     "( 0  0   0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector multiplication assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec( 4UL, 0 );
      vec[0] =  2;
      vec[1] = -4;

      col2 *= vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector multiplication assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 *= vec;

      checkSize    ( col2 , 4UL );
      checkNonZeros( col2 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major scalar multiplication assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      col2 *= 3;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != -6 || col2[1] != 0 || col2[2] != -9 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -9 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseColumn division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the SparseColumn class
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

      CT col2 = column( mat_, 2UL );

      col2 /= 0.5;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != -6 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -6 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar division assignment
   //=====================================================================================

   {
      test_ = "Column-major scalar division assignment";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      col2 /= 0.5;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != -6 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -6 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseColumn subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the SparseColumn class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSubscript()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::operator[]";

      initialize();

      CT col2 = column( mat_, 2UL );

      // Writing the first element
      col2[1] = 9;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the second element
      col2[2] = 0;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 2UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the third element
      col2[3] = -8;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != -8 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::operator[]";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      // Writing the first element
      col2[1] = 9;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the second element
      col2[2] = 0;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 2UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the third element
      col2[3] = -8;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != -8 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseColumn iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the SparseColumn class template.
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

         CT col0 = column( mat_, 0UL );
         const size_t number( col0.end() - col0.begin() );

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

         CT col1 = column( mat_, 1UL );
         const size_t number( col1.end() - col1.begin() );

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

         CT col2 = column( mat_, 2UL );
         const size_t number( col2.end() - col2.begin() );

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

         CT col3 = column( mat_, 3UL );
         const size_t number( col3.end() - col3.begin() );

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

         CT col4 = column( mat_, 4UL );
         const size_t number( col4.end() - col4.begin() );

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

         CT col2 = column( mat_, 2UL );
         CT::ConstIterator it ( col2.cbegin() );
         CT::ConstIterator end( col2.cend() );

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

         CT col4 = column( mat_, 4UL );
         int value = 6;

         for( CT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it = value++;
         }

         if( col4[0] != 6 || col4[1] != 7 || col4[2] != 8 || col4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) != 6 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 7 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 8 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5  8 )\n"
                                        "( 0  0  0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         CT col4 = column( mat_, 4UL );
         int value = 2;

         for( CT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it += value++;
         }

         if( col4[0] != 8 || col4[1] != 10 || col4[2] != 12 || col4[3] != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  8 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 10 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 12 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  8 )\n"
                                        "( 0  1  0  4 10 )\n"
                                        "( 0  0 -3  5 12 )\n"
                                        "( 0  0  0 -6 14 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         CT col4 = column( mat_, 4UL );
         int value = 2;

         for( CT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it -= value++;
         }

         if( col4[0] != 6 || col4[1] != 7 || col4[2] != 8 || col4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) != 6 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 7 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 8 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5  8 )\n"
                                        "( 0  0  0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         CT col4 = column( mat_, 4UL );
         int value = 1;

         for( CT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it *= value++;
         }

         if( col4[0] != 6 || col4[1] != 14 || col4[2] != 24 || col4[3] != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  6 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 14 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 24 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4 14 )\n"
                                        "( 0  0 -3  5 24 )\n"
                                        "( 0  0  0 -6 36 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         CT col4 = column( mat_, 4UL );

         for( CT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it /= 2;
         }

         if( col4[0] != 3 || col4[1] != 7 || col4[2] != 12 || col4[3] != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 3 7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  3 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) !=  7 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 12 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  3 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5 12 )\n"
                                        "( 0  0  0 -6 18 )\n";
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

         TCT col0 = column( tmat_, 0UL );
         const size_t number( col0.end() - col0.begin() );

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

         TCT col1 = column( tmat_, 1UL );
         const size_t number( col1.end() - col1.begin() );

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

         TCT col2 = column( tmat_, 2UL );
         const size_t number( col2.end() - col2.begin() );

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

         TCT col3 = column( tmat_, 3UL );
         const size_t number( col3.end() - col3.begin() );

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

         TCT col4 = column( tmat_, 4UL );
         const size_t number( col4.end() - col4.begin() );

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

         TCT col2 = column( tmat_, 2UL );
         TCT::ConstIterator it ( col2.cbegin() );
         TCT::ConstIterator end( col2.cend() );

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

         TCT col4 = column( tmat_, 4UL );
         int value = 6;

         for( TCT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it = value++;
         }

         if( col4[0] != 6 || col4[1] != 7 || col4[2] != 8 || col4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) != 6 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 7 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 8 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5  8 )\n"
                                        "( 0  0  0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         TCT col4 = column( tmat_, 4UL );
         int value = 2;

         for( TCT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it += value++;
         }

         if( col4[0] != 8 || col4[1] != 10 || col4[2] != 12 || col4[3] != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  8 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 10 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 12 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  8 )\n"
                                        "( 0  1  0  4 10 )\n"
                                        "( 0  0 -3  5 12 )\n"
                                        "( 0  0  0 -6 14 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         TCT col4 = column( tmat_, 4UL );
         int value = 2;

         for( TCT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it -= value++;
         }

         if( col4[0] != 6 || col4[1] != 7 || col4[2] != 8 || col4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) != 6 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 7 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 8 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5  8 )\n"
                                        "( 0  0  0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         TCT col4 = column( tmat_, 4UL );
         int value = 1;

         for( TCT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it *= value++;
         }

         if( col4[0] != 6 || col4[1] != 14 || col4[2] != 24 || col4[3] != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  6 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 14 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 24 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4 14 )\n"
                                        "( 0  0 -3  5 24 )\n"
                                        "( 0  0  0 -6 36 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         TCT col4 = column( tmat_, 4UL );

         for( TCT::Iterator it=col4.begin(); it!=col4.end(); ++it ) {
            *it /= 2;
         }

         if( col4[0] != 3 || col4[1] != 7 || col4[2] != 12 || col4[3] != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 3 7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  3 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) !=  7 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 12 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  3 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5 12 )\n"
                                        "( 0  0  0 -6 18 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the nonZeros member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of SparseColumn. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::nonZeros()";

      initialize();

      // Initialization check
      CT col3 = column( mat_, 3UL );

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 3UL );

      if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 0 4 5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse column
      col3[2] = 0;

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 2UL );

      if( col3[0] != 0 || col3[1] != 4 || col3[2] != 0 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 0 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      mat_(0,3) = 5;

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 3UL );

      if( col3[0] != 5 || col3[1] != 4 || col3[2] != 0 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 5 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::nonZeros()";

      initialize();

      // Initialization check
      TCT col3 = column( tmat_, 3UL );

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 3UL );

      if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 0 4 5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse column
      col3[2] = 0;

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 2UL );

      if( col3[0] != 0 || col3[1] != 4 || col3[2] != 0 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 0 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      tmat_(0,3) = 5;

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 3UL );

      if( col3[0] != 5 || col3[1] != 4 || col3[2] != 0 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 5 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reset member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of SparseColumn. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::reset()";

      initialize();

      // Resetting the 0th column
      {
         CT col0 = column( mat_, 0UL );
         col0.reset();

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  0UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 10UL );

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st column
      {
         CT col1 = column( mat_, 1UL );
         col1.reset();

         checkSize    ( col1, 4UL );
         checkNonZeros( col1, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 9UL );

         if( col1[0] != 0 || col1[1] != 0 || col1[2] != 0 || col1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st column failed\n"
                << " Details:\n"
                << "   Result:\n" << col1 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd column
      {
         CT col2 = column( mat_, 2UL );
         col2.reset();

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 7UL );

         if( col2[0] != 0 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 3rd column
      {
         CT col3 = column( mat_, 3UL );
         col3.reset();

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 4UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 0 || col3[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 4th column
      {
         CT col4 = column( mat_, 4UL );
         col4.reset();

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 0UL );

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 4th column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::reset()";

      initialize();

      // Resetting the 0th column
      {
         TCT col0 = column( tmat_, 0UL );
         col0.reset();

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  0UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 10UL );

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st column
      {
         TCT col1 = column( tmat_, 1UL );
         col1.reset();

         checkSize    ( col1 , 4UL );
         checkNonZeros( col1 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 9UL );

         if( col1[0] != 0 || col1[1] != 0 || col1[2] != 0 || col1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st column failed\n"
                << " Details:\n"
                << "   Result:\n" << col1 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd column
      {
         TCT col2 = column( tmat_, 2UL );
         col2.reset();

         checkSize    ( col2 , 4UL );
         checkNonZeros( col2 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( col2[0] != 0 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 3rd column
      {
         TCT col3 = column( tmat_, 3UL );
         col3.reset();

         checkSize    ( col3 , 4UL );
         checkNonZeros( col3 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 4UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 0 || col3[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 4th column
      {
         TCT col4 = column( tmat_, 4UL );
         col4.reset();

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 0UL );

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 4th column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the append member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the append member function of SparseColumn. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::append()";

      MT mat( 9UL, 3UL );

      CT col1 = column( mat, 1UL );
      col1.reserve( 4UL );

      // Appending one non-zero element
      col1.append( 1UL, 1 );

      checkSize    ( col1, 9UL );
      checkCapacity( col1, 4UL );
      checkNonZeros( col1, 1UL );

      if( col1[1] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Appending three more non-zero elements
      col1.append( 3UL, 2 );
      col1.append( 4UL, 3 );
      col1.append( 8UL, 4 );

      checkSize    ( col1, 9UL );
      checkCapacity( col1, 4UL );
      checkNonZeros( col1, 4UL );

      if( col1[1] != 1 || col1[3] != 2 || col1[4] != 3 || col1[8] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 0 2 3 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::append()";

      TMT mat( 9UL, 3UL );

      TCT col1 = column( mat, 1UL );
      col1.reserve( 4UL );

      // Appending one non-zero element
      col1.append( 1UL, 1 );

      checkSize    ( col1, 9UL );
      checkCapacity( col1, 4UL );
      checkNonZeros( col1, 1UL );

      if( col1[1] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Appending three more non-zero elements
      col1.append( 3UL, 2 );
      col1.append( 4UL, 3 );
      col1.append( 8UL, 4 );

      checkSize    ( col1, 9UL );
      checkCapacity( col1, 4UL );
      checkNonZeros( col1, 4UL );

      if( col1[1] != 1 || col1[3] != 2 || col1[4] != 3 || col1[8] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 0 2 3 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the insert member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the insert member function of SparseColumn. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::insert()";

      initialize();

      CT col0 = column( mat_, 0UL );

      // Inserting a non-zero element at the end of the column
      {
         CT::Iterator pos = col0.insert( 3UL, 1 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  1UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
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

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the beginning of the column
      {
         CT::Iterator pos = col0.insert( 0UL, 2 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  2UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the center of the column
      {
         CT::Iterator pos = col0.insert( 2UL, 3 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  3UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 3 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         col0.insert( 3UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 2 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::insert()";

      initialize();

      TCT col0 = column( tmat_, 0UL );

      // Inserting a non-zero element at the end of the column
      {
         TCT::Iterator pos = col0.insert( 3UL, 1 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
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

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the beginning of the column
      {
         TCT::Iterator pos = col0.insert( 0UL, 2 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  2UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the center of the column
      {
         TCT::Iterator pos = col0.insert( 2UL, 3 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  3UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 3 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         col0.insert( 3UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 2 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the erase member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the erase member function of SparseColumn. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::erase( size_t )";

      initialize();

      CT col4 = column( mat_, 4UL );

      // Erasing the non-zero element at the end of the column
      col4.erase( 3UL );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 3UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 7 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the column
      col4.erase( 0UL );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 2UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 8UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the column
      col4.erase( 2UL );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 7UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      col4.erase( 3UL );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 7UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::erase( Iterator )";

      initialize();

      CT col4 = column( mat_, 4UL );

      // Erasing the non-zero element at the end of the column
      {
         CT::Iterator pos = col4.erase( col4.find( 3UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 3UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 9UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 7 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the column
      {
         CT::Iterator pos = col4.erase( col4.find( 0UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 2UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
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

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the column
      {
         CT::Iterator pos = col4.erase( col4.find( 2UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 1UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 7UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         CT::Iterator pos = col4.erase( col4.find( 3UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 1UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 7UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::erase( Iterator, Iterator )";

      initialize();

      // Erasing the 2nd column
      {
         CT col2 = column( mat_, 2UL );

         CT::Iterator pos = col2.erase( col2.begin(), col2.end() );

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 8UL );

         if( pos != col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col2[0] != 0 || col2[1] != -0 || col2[2] != 0 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 4th column
      {
         CT col4 = column( mat_, 4UL );

         CT::Iterator pos = col4.erase( col4.begin(), col4.find( 2UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 2UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
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

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 9 || col4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 4th column
      {
         CT col4 = column( mat_, 4UL );

         CT::Iterator pos = col4.erase( col4.find( 2UL ), col4.end() );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 4UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         CT col3 = column( mat_, 3UL );

         CT::Iterator pos = col3.erase( col3.find( 1UL ), col3.find( 1UL ) );

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 3UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 4UL );

         if( pos != col3.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::erase( size_t )";

      initialize();

      TCT col4 = column( tmat_, 4UL );

      // Erasing the non-zero element at the end of the column
      col4.erase( 3UL );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 3UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 7 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the column
      col4.erase( size_t(0) );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 8UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the column
      col4.erase( 2UL );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      col4.erase( 3UL );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::erase( size_t )";

      initialize();

      TCT col4 = column( tmat_, 4UL );

      // Erasing the non-zero element at the end of the column
      {
         TCT::Iterator pos = col4.erase( col4.find( 3UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 3UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 9UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 7 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the column
      {
         TCT::Iterator pos = col4.erase( col4.find( 0UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 2UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
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

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the column
      {
         TCT::Iterator pos = col4.erase( col4.find( 2UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 1UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         TCT::Iterator pos = col4.erase( col4.find( 3UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 1UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::erase( Iterator, Iterator )";

      initialize();

      // Erasing the 2nd column
      {
         TCT col2 = column( tmat_, 2UL );

         TCT::Iterator pos = col2.erase( col2.begin(), col2.end() );

         checkSize    ( col2 , 4UL );
         checkNonZeros( col2 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 8UL );

         if( pos != col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col2[0] != 0 || col2[1] != -0 || col2[2] != 0 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 4th column
      {
         TCT col4 = column( tmat_, 4UL );

         TCT::Iterator pos = col4.erase( col4.begin(), col4.find( 2UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 2UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
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

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 9 || col4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 4th column
      {
         TCT col4 = column( tmat_, 4UL );

         TCT::Iterator pos = col4.erase( col4.find( 2UL ), col4.end() );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 4UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         TCT col3 = column( tmat_, 3UL );

         TCT::Iterator pos = col3.erase( col3.find( 1UL ), col3.find( 1UL ) );

         checkSize    ( col3 , 4UL );
         checkNonZeros( col3 , 3UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 4UL );

         if( pos != col3.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reserve member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reserve member function of SparseColumn. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::reserve()";

      MT mat( 20UL, 3UL );

      CT col0 = column( mat, 0UL );

      // Increasing the capacity of the column
      col0.reserve( 10UL );

      checkSize    ( col0, 20UL );
      checkCapacity( col0, 10UL );
      checkNonZeros( col0,  0UL );

      // Further increasing the capacity of the column
      col0.reserve( 15UL );

      checkSize    ( col0, 20UL );
      checkCapacity( col0, 15UL );
      checkNonZeros( col0,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::reserve()";

      TMT mat( 20UL, 3UL );

      TCT col0 = column( mat, 0UL );

      // Increasing the capacity of the column
      col0.reserve( 10UL );

      checkSize    ( col0, 20UL );
      checkCapacity( col0, 10UL );
      checkNonZeros( col0,  0UL );

      // Further increasing the capacity of the column
      col0.reserve( 15UL );

      checkSize    ( col0, 20UL );
      checkCapacity( col0, 15UL );
      checkNonZeros( col0,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of SparseColumn. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScale()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::scale()";

      initialize();

      // Integral scaling the 3rd column
      {
         CT col3 = column( mat_, 3UL );
         col3.scale( 3 );

         checkSize    ( col3,  4UL );
         checkNonZeros( col3,  3UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 10UL );

         if( col3[0] != 0 || col3[1] != 12 || col3[2] != 15 || col3[3] != -18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 12 15 -18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=   0 || mat_(0,4) !=  7 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  12 || mat_(1,4) != -8 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  15 || mat_(2,4) !=  9 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -18 || mat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2   0  7 )\n"
                                        "( 0  1  0  12 -8 )\n"
                                        "( 0  0 -3  15  9 )\n"
                                        "( 0  0  0 -18 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Floating point scaling the 3rd column
      {
         CT col3 = column( mat_, 3UL );
         col3.scale( 0.5 );

         checkSize    ( col3,  4UL );
         checkNonZeros( col3,  3UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 10UL );

         if( col3[0] != 0 || col3[1] != 6 || col3[2] != 7 || col3[3] != -9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 6 7 -9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  6 || mat_(1,4) != -8 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  7 || mat_(2,4) !=  9 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -9 || mat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  1  0  6 -8 )\n"
                                        "( 0  0 -3  7  9 )\n"
                                        "( 0  0  0 -9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::scale()";

      initialize();

      // Integral scaling the 3rd column
      {
         TCT col3 = column( tmat_, 3UL );
         col3.scale( 3 );

         checkSize    ( col3 ,  4UL );
         checkNonZeros( col3 ,  3UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 10UL );

         if( col3[0] != 0 || col3[1] != 12 || col3[2] != 15 || col3[3] != -18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 12 15 -18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=   0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  12 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  15 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -18 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2   0  7 )\n"
                                        "( 0  1  0  12 -8 )\n"
                                        "( 0  0 -3  15  9 )\n"
                                        "( 0  0  0 -18 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Floating point scaling the 3rd column
      {
         TCT col3 = column( tmat_, 3UL );
         col3.scale( 0.5 );

         checkSize    ( col3 ,  4UL );
         checkNonZeros( col3 ,  3UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 10UL );

         if( col3[0] != 0 || col3[1] != 6 || col3[2] != 7 || col3[3] != -9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 6 7 -9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  6 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  7 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -9 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  1  0  6 -8 )\n"
                                        "( 0  0 -3  7  9 )\n"
                                        "( 0  0  0 -9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the find member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the find member function of SparseColumn. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::find()";

      initialize();

      CT col2 = column( mat_, 2UL );

      // Searching for the first element
      {
         CT::Iterator pos = col2.find( 0UL );

         if( pos == col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col2 << "\n";
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
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         CT::Iterator pos = col2.find( 2UL );

         if( pos == col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col2 << "\n";
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
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         CT::Iterator pos = col2.find( 1UL );

         if( pos != col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::find()";

      initialize();

      TCT col2 = column( tmat_, 2UL );

      // Searching for the first element
      {
         TCT::Iterator pos = col2.find( 0UL );

         if( pos == col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col2 << "\n";
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
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         TCT::Iterator pos = col2.find( 2UL );

         if( pos == col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col2 << "\n";
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
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         TCT::Iterator pos = col2.find( 1UL );

         if( pos != col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the lowerBound member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the lowerBound member function of SparseColumn. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::lowerBound()";

      initialize();

      CT col1 = column( mat_, 1UL );

      // Determining the lower bound for index 0
      {
         CT::Iterator pos = col1.lowerBound( 0UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 1
      {
         CT::Iterator pos = col1.lowerBound( 1UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 2
      {
         CT::Iterator pos = col1.lowerBound( 2UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::lowerBound()";

      initialize();

      TCT col1 = column( tmat_, 1UL );

      // Determining the lower bound for index 0
      {
         TCT::Iterator pos = col1.lowerBound( 0UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 1
      {
         TCT::Iterator pos = col1.lowerBound( 1UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 2
      {
         TCT::Iterator pos = col1.lowerBound( 2UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the upperBound member function of SparseColumn.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the upperBound member function of SparseColumn. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseColumn::upperBound()";

      initialize();

      CT col1 = column( mat_, 1UL );

      // Determining the upper bound for index 0
      {
         CT::Iterator pos = col1.upperBound( 0UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 1
      {
         CT::Iterator pos = col1.upperBound( 1UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 2
      {
         CT::Iterator pos = col1.upperBound( 2UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseColumn::upperBound()";

      initialize();

      TCT col1 = column( tmat_, 1UL );

      // Determining the upper bound for index 0
      {
         TCT::Iterator pos = col1.upperBound( 0UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 1
      {
         TCT::Iterator pos = col1.upperBound( 1UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 2
      {
         TCT::Iterator pos = col1.upperBound( 2UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDefault function with the SparseColumn class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDefault function with the SparseColumn class template.
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

      // isDefault with default column
      {
         CT col0 = column( mat_, 0UL );

         if( isDefault( col0 ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default column
      {
         CT col1 = column( mat_, 1UL );

         if( isDefault( col1 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col1 << "\n";
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

      // isDefault with default column
      {
         TCT col0 = column( tmat_, 0UL );

         if( isDefault( col0 ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default column
      {
         TCT col1 = column( tmat_, 1UL );

         if( isDefault( col1 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isnan function with the SparseColumn class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isnan function with the SparseColumn class template.
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
      typedef blaze::SparseColumn<MatrixType>                 ColumnType;

      MatrixType mat( mat_ );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  5UL );
      checkNonZeros( mat, 10UL );

      // isnan with empty column
      {
         ColumnType col0 = column( mat, 0UL );

         checkSize    ( col0, 4UL );
         checkNonZeros( col0, 0UL );

         if( blaze::isnan( col0 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with partially filled column
      {
         ColumnType col2 = column( mat, 2UL );

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 2UL );

         if( blaze::isnan( col2 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with fully filled column
      {
         ColumnType col4 = column( mat, 4UL );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 4UL );

         if( blaze::isnan( col4 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col4 << "\n";
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
      typedef blaze::SparseColumn<MatrixType>                    ColumnType;

      MatrixType mat( mat_ );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  5UL );
      checkNonZeros( mat, 10UL );

      // isnan with empty column
      {
         ColumnType col0 = column( mat, 0UL );

         checkSize    ( col0, 4UL );
         checkNonZeros( col0, 0UL );

         if( blaze::isnan( col0 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with partially filled column
      {
         ColumnType col2 = column( mat, 2UL );

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 2UL );

         if( blaze::isnan( col2 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with fully filled column
      {
         ColumnType col4 = column( mat, 4UL );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 4UL );

         if( blaze::isnan( col4 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col4 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the min function with the SparseColumn class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the min function used with the SparseColumn class
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

      // Computing the minimum of the 0th column
      {
         const int minimum = min( column( mat_, 0UL ) );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 0th column failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 1st column
      {
         const int minimum = min( column( mat_, 1UL ) );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 1st column failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 2nd column
      {
         const int minimum = min( column( mat_, 2UL ) );

         if( minimum != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 2nd column failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 3rd column
      {
         const int minimum = min( column( mat_, 3UL ) );

         if( minimum != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 3rd column failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 4th column
      {
         const int minimum = min( column( mat_, 4UL ) );

         if( minimum != -8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 4th column failed\n"
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

      // Computing the minimum of the 0th column
      {
         const int minimum = min( column( tmat_, 0UL ) );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 0th column failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 1st column
      {
         const int minimum = min( column( tmat_, 1UL ) );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 1st column failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 2nd column
      {
         const int minimum = min( column( tmat_, 2UL ) );

         if( minimum != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 2nd column failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 3rd column
      {
         const int minimum = min( column( tmat_, 3UL ) );

         if( minimum != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 3rd column failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the minimum of the 4th column
      {
         const int minimum = min( column( tmat_, 4UL ) );

         if( minimum != -8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Minimum computation for 4th column failed\n"
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
/*!\brief Test of the max function with the SparseColumn class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the max function used with the SparseColumn class
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

      // Computing the maximum of the 0th column
      {
         const int maximum = max( column( mat_, 0UL ) );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 0th column failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 1st column
      {
         const int maximum = max( column( mat_, 1UL ) );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 1st column failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 2nd column
      {
         const int maximum = max( column( mat_, 2UL ) );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 2nd column failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 3rd column
      {
         const int maximum = max( column( mat_, 3UL ) );

         if( maximum != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 3rd column failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 4th column
      {
         const int maximum = max( column( mat_, 4UL ) );

         if( maximum != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 4th column failed\n"
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

      // Computing the maximum of the 0th column
      {
         const int maximum = max( column( tmat_, 0UL ) );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 0th column failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 1st column
      {
         const int maximum = max( column( tmat_, 1UL ) );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 1st column failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 2nd column
      {
         const int maximum = max( column( tmat_, 2UL ) );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 2nd column failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 3rd column
      {
         const int maximum = max( column( tmat_, 3UL ) );

         if( maximum != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 3rd column failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Computing the maximum of the 4th column
      {
         const int maximum = max( column( tmat_, 4UL ) );

         if( maximum != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Maximum computation for 4th column failed\n"
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
/*!\brief Test of the subvector function with the SparseColumn class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subvector function used with the SparseColumn class
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

      typedef blaze::SparseSubvector<CT>  SubvectorType;

      CT col1 = column( mat_, 1UL );
      SubvectorType sv = subvector( col1, 0UL, 4UL );

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

      typedef blaze::SparseSubvector<TCT>  SubvectorType;

      TCT col1 = column( tmat_, 1UL );
      SubvectorType sv = subvector( col1, 0UL, 4UL );

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
   mat_(0,2) = -2;
   mat_(2,2) = -3;
   mat_(1,3) =  4;
   mat_(2,3) =  5;
   mat_(3,3) = -6;
   mat_(0,4) =  7;
   mat_(1,4) = -8;
   mat_(2,4) =  9;
   mat_(3,4) = 10;

   // Initializing the column-major compressed matrix
   tmat_.reset();
   tmat_(1,1) =  1;
   tmat_(0,2) = -2;
   tmat_(2,2) = -3;
   tmat_(1,3) =  4;
   tmat_(2,3) =  5;
   tmat_(3,3) = -6;
   tmat_(0,4) =  7;
   tmat_(1,4) = -8;
   tmat_(2,4) =  9;
   tmat_(3,4) = 10;
}
//*************************************************************************************************

} // namespace sparsecolumn

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
   std::cout << "   Running SparseColumn class test..." << std::endl;

   try
   {
      RUN_SPARSECOLUMN_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SparseColumn class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
