//=================================================================================================
/*!
//  \file src/mathtest/submatrix/DenseUnalignedTest.cpp
//  \brief Source file for the Submatrix dense unaligned test
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
#include <memory>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/Views.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blazetest/mathtest/submatrix/DenseUnalignedTest.h>


namespace blazetest {

namespace mathtest {

namespace submatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Submatrix dense unaligned test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseUnalignedTest::DenseUnalignedTest()
   : mat_ ( 5UL, 4UL )
   , tmat_( 4UL, 5UL )
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
   testTranspose();
   testCTranspose();
   testIsDefault();
   testIsSame();
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
/*!\brief Test of the Submatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Submatrix specialization. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testConstructors()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major Submatrix constructor";

      initialize();

      for( size_t row=0UL; row<mat_.rows(); ++row ) {
         for( size_t column=0UL; column<mat_.columns(); ++column ) {
            for( size_t m=0UL; (row+m)<mat_.rows(); ++m ) {
               for( size_t n=0UL; (column+n)<mat_.columns(); ++n )
               {
                  SMT sm = blaze::submatrix( mat_, row, column, m, n );

                  for( size_t i=0UL; i<m; ++i ) {
                     for( size_t j=0UL; j<n; ++j )
                     {
                        if( sm(i,j) != mat_(row+i,column+j) ) {
                           std::ostringstream oss;
                           oss << " Test: " << test_ << "\n"
                               << " Error: Setup of dense submatrix failed\n"
                               << " Details:\n"
                               << "   Index of first row    = " << row << "\n"
                               << "   Index of first column = " << column << "\n"
                               << "   Number of rows        = " << m << "\n"
                               << "   Number of columns     = " << n << "\n"
                               << "   Submatrix:\n" << sm << "\n"
                               << "   Matrix:\n" << mat_ << "\n";
                           throw std::runtime_error( oss.str() );
                        }
                     }
                  }
               }
            }
         }
      }

      try {
         SMT sm = blaze::submatrix( mat_, 2UL, 2UL, 4UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm = blaze::submatrix( mat_, 2UL, 2UL, 2UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm = blaze::submatrix( mat_, 5UL, 2UL, 2UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm = blaze::submatrix( mat_, 2UL, 4UL, 2UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major Submatrix constructor";

      initialize();

      for( size_t column=0UL; column<tmat_.columns(); ++column ) {
         for( size_t row=0UL; row<tmat_.rows(); ++row ) {
            for( size_t n=0UL; (column+n)<tmat_.columns(); ++n ) {
               for( size_t m=0UL; (row+m)<tmat_.rows(); ++m )
               {
                  OSMT sm = blaze::submatrix( tmat_, row, column, m, n );

                  for( size_t j=0UL; j<n; ++j ) {
                     for( size_t i=0UL; i<m; ++i )
                     {
                        if( sm(i,j) != tmat_(row+i,column+j) ) {
                           std::ostringstream oss;
                           oss << " Test: " << test_ << "\n"
                               << " Error: Setup of dense submatrix failed\n"
                               << " Details:\n"
                               << "   Index of first row    = " << row << "\n"
                               << "   Index of first column = " << column << "\n"
                               << "   Number of rows        = " << m << "\n"
                               << "   Number of columns     = " << n << "\n"
                               << "   Submatrix:\n" << sm << "\n"
                               << "   Matrix:\n" << tmat_ << "\n";
                           throw std::runtime_error( oss.str() );
                        }
                     }
                  }
               }
            }
         }
      }

      try {
         OSMT sm = blaze::submatrix( tmat_, 2UL, 2UL, 3UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         OSMT sm = blaze::submatrix( tmat_, 2UL, 2UL, 2UL, 4UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         OSMT sm = blaze::submatrix( tmat_, 4UL, 2UL, 2UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         OSMT sm = blaze::submatrix( tmat_, 2UL, 5UL, 2UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testAssignment()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix homogeneous assignment";

      initialize();

      // Assigning to a 2x3 submatrix
      {
         SMT sm = blaze::submatrix( mat_, 0UL, 1UL, 2UL, 3UL );
         sm = 12;

         checkRows    ( sm  ,  2UL );
         checkColumns ( sm  ,  3UL );
         checkNonZeros( sm  ,  6UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 15UL );

         if( sm(0,0) != 12 || sm(0,1) != 12 || sm(0,2) != 12 ||
             sm(1,0) != 12 || sm(1,1) != 12 || sm(1,2) != 12 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 12 12 )\n( 12 12 12 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) != 12 || mat_(0,2) != 12 || mat_(0,3) != 12 ||
             mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) != 12 || mat_(1,3) != 12 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0 12 12 12 )\n"
                                        "(  0 12 12 12 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  4  5 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Assigning to a 3x2 submatrix
      {
         SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 2UL );
         sm = 15;

         checkRows    ( sm  ,  3UL );
         checkColumns ( sm  ,  2UL );
         checkNonZeros( sm  ,  6UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 18UL );

         if( sm(0,0) != 15 || sm(1,1) != 15 ||
             sm(1,0) != 15 || sm(1,1) != 15 ||
             sm(2,0) != 15 || sm(2,1) != 15 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 15 15 )\n( 15 15 )\n( 15 15 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) != 12 || mat_(0,2) != 12 || mat_(0,3) != 12 ||
             mat_(1,0) != 15 || mat_(1,1) != 15 || mat_(1,2) != 12 || mat_(1,3) != 12 ||
             mat_(2,0) != 15 || mat_(2,1) != 15 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) != 15 || mat_(3,1) != 15 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0 12 12 12 )\n"
                                        "( 15 15 12 12 )\n"
                                        "( 15 15 -3  0 )\n"
                                        "( 15 15  5 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major initializer list assignment (complete list)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );
      sm = { { 1, 2, 3 }, { 4, 5, 6 } };

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 13UL );

      if( sm(0,0) != 1 || sm(0,1) != 2 || sm(0,2) != 3 ||
          sm(1,0) != 4 || sm(1,1) != 5 || sm(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 1 || mat_(1,1) !=  2 || mat_(1,2) != 3 || mat_(1,3) !=  0 ||
          mat_(2,0) != 4 || mat_(2,1) !=  5 || mat_(2,2) != 6 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  4 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 1  2  3  0 )\n"
                                     "( 4  5  6  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major initializer list assignment (incomplete list)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );
      sm = { { 1 }, { 4, 5, 6 } };

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) != 1 || sm(0,1) != 0 || sm(0,2) != 0 ||
          sm(1,0) != 4 || sm(1,1) != 5 || sm(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 1 || mat_(1,1) !=  0 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 4 || mat_(2,1) !=  5 || mat_(2,2) != 6 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  4 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 1  0  0  0 )\n"
                                     "( 4  5  6  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix copy assignment (no aliasing)";

      initialize();

      MT mat( 5UL, 4UL, 0 );
      mat(1,0) = 11;
      mat(2,0) = 12;
      mat(2,2) = 13;

      SMT sm = blaze::submatrix( mat, 1UL, 0UL, 2UL, 3UL );
      sm = blaze::submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );
      checkRows    ( mat ,  5UL );
      checkColumns ( mat ,  4UL );
      checkNonZeros( mat ,  4UL );

      if( sm(0,0) != 0 || sm(0,1) != -3 || sm(0,2) !=  0 ||
          sm(1,0) != 4 || sm(1,1) !=  5 || sm(1,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -3  0 )\n( 4  5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  0 || mat(0,2) !=  0 || mat(0,3) != 0 ||
          mat(1,0) != 0 || mat(1,1) != -3 || mat(1,2) !=  0 || mat(1,3) != 0 ||
          mat(2,0) != 4 || mat(2,1) !=  5 || mat(2,2) != -6 || mat(2,3) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) !=  0 || mat(3,3) != 0 ||
          mat(4,0) != 0 || mat(4,1) !=  0 || mat(4,2) !=  0 || mat(4,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0 -3  0  0 )\n"
                                     "( 4  5 -6  0 )\n"
                                     "( 0  0  0  0 )\n"
                                     "( 0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix copy assignment (aliasing)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );
      sm = blaze::submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) != 0 || sm(0,1) != -3 || sm(0,2) !=  0 ||
          sm(1,0) != 4 || sm(1,1) !=  5 || sm(1,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -3  0 )\n( 4  5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -3 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 4 || mat_(2,1) !=  5 || mat_(2,2) != -6 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0 -3  0  0 )\n"
                                     "( 4  5 -6  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix assignment (mixed type)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {  0, 11,  0 },
                                                      { 12, 13, 14 } };

      sm = mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 11 || sm(0,2) !=  0 ||
          sm(1,0) != 12 || sm(1,1) != 13 || sm(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 11  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 11 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 12 || mat_(2,1) != 13 || mat_(2,2) != 14 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 11  0  0 )\n"
                                     "( 12 13 14  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix assignment (aligned/padded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm = mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 11 || sm(0,2) !=  0 ||
          sm(1,0) != 12 || sm(1,1) != 13 || sm(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 11  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 11 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 12 || mat_(2,1) != 13 || mat_(2,2) != 14 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 11  0  0 )\n"
                                     "( 12 13 14  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 3UL );
      mat = 0;
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm = mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 11 || sm(0,2) !=  0 ||
          sm(1,0) != 12 || sm(1,1) != 13 || sm(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 11  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 11 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 12 || mat_(2,1) != 13 || mat_(2,2) != 14 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 11  0  0 )\n"
                                     "( 12 13 14  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (mixed type)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {  0, 11,  0 },
                                                         { 12, 13, 14 } };

      sm = mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 11 || sm(0,2) !=  0 ||
          sm(1,0) != 12 || sm(1,1) != 13 || sm(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 11  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 11 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 12 || mat_(2,1) != 13 || mat_(2,2) != 14 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 11  0  0 )\n"
                                     "( 12 13 14  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (aligned/padded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm = mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 11 || sm(0,2) !=  0 ||
          sm(1,0) != 12 || sm(1,1) != 13 || sm(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 11  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 11 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 12 || mat_(2,1) != 13 || mat_(2,2) != 14 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 11  0  0 )\n"
                                     "( 12 13 14  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 3UL );
      mat = 0;
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm = mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 11 || sm(0,2) !=  0 ||
          sm(1,0) != 12 || sm(1,1) != 13 || sm(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 11  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 11 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 12 || mat_(2,1) != 13 || mat_(2,2) != 14 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 11  0  0 )\n"
                                     "( 12 13 14  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix assignment";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm = mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 11 || sm(0,2) !=  0 ||
          sm(1,0) != 12 || sm(1,1) != 13 || sm(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 11  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 11 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 12 || mat_(2,1) != 13 || mat_(2,2) != 14 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 11  0  0 )\n"
                                     "( 12 13 14  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix assignment";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm = mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 11 || sm(0,2) !=  0 ||
          sm(1,0) != 12 || sm(1,1) != 13 || sm(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 11  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 11 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 12 || mat_(2,1) != 13 || mat_(2,2) != 14 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 11  0  0 )\n"
                                     "( 12 13 14  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix homogeneous assignment";

      initialize();

      // Assigning to a 3x2 submatrix
      {
         OSMT sm = blaze::submatrix( tmat_, 1UL, 0UL, 3UL, 2UL );
         sm = 12;

         checkRows    ( sm   ,  3UL );
         checkColumns ( sm   ,  2UL );
         checkNonZeros( sm   ,  6UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 15UL );

         if( sm(0,0) != 12 || sm(0,1) != 12 ||
             sm(1,0) != 12 || sm(1,1) != 12 ||
             sm(2,0) != 12 || sm(2,1) != 12 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 12 )\n( 12 12 )\n( 12 12 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 12 || tmat_(1,1) != 12 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 12 || tmat_(2,1) != 12 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 12 || tmat_(3,1) != 12 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0  0 -2  0  7 )\n"
                                        "( 12 12  0  4 -8 )\n"
                                        "( 12 12 -3  5  9 )\n"
                                        "( 12 12  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Assigning to a 2x3 submatrix
      {
         OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 3UL );
         sm = 15;

         checkRows    ( sm   ,  2UL );
         checkColumns ( sm   ,  3UL );
         checkNonZeros( sm   ,  6UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 18UL );

         if( sm(0,0) != 15 || sm(0,1) != 15 || sm(0,2) != 15 ||
             sm(1,0) != 15 || sm(1,1) != 15 || sm(1,2) != 15 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 15 15 15 )\n( 15 15 15 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) !=  0 || tmat_(0,1) != 15 || tmat_(0,2) != 15 || tmat_(0,3) != 15 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 12 || tmat_(1,1) != 15 || tmat_(1,2) != 15 || tmat_(1,3) != 15 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 12 || tmat_(2,1) != 12 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 12 || tmat_(3,1) != 12 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n(  0 15 15 15  7 )\n"
                                        "( 12 15 15 15 -8 )\n"
                                        "( 12 12 -3  5  9 )\n"
                                        "( 12 12  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major initializer list assignment (complete list)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );
      sm = { { 1, 2 }, { 3, 4 }, { 5, 6 } };

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 13UL );

      if( sm(0,0) != 1 || sm(0,1) != 2 ||
          sm(1,0) != 3 || sm(1,1) != 4 ||
          sm(2,0) != 5 || sm(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 2 )\n( 3 4 )\n( 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 1 || tmat_(0,2) != 2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 3 || tmat_(1,2) != 4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 5 || tmat_(2,2) != 6 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  1  2  0  7 )\n"
                                     "( 0  3  4  4 -8 )\n"
                                     "( 0  5  6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major initializer list assignment (incomplete list)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );
      sm = { { 1 }, { 3 }, { 5, 6 } };

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) != 1 || sm(0,1) != 0 ||
          sm(1,0) != 3 || sm(1,1) != 0 ||
          sm(2,0) != 5 || sm(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 )\n( 3 0 )\n( 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 1 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 3 || tmat_(1,2) != 0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 5 || tmat_(2,2) != 6 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  1  0  0  7 )\n"
                                     "( 0  3  0  4 -8 )\n"
                                     "( 0  5  6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix copy assignment (no aliasing)";

      initialize();

      OMT mat( 4UL, 5UL, 0 );
      mat(0,1) = 11;
      mat(0,2) = 12;
      mat(2,2) = 13;

      OSMT sm = blaze::submatrix( mat, 0UL, 1UL, 3UL, 2UL );
      sm = blaze::submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );
      checkRows    ( mat  ,  4UL );
      checkColumns ( mat  ,  5UL );
      checkNonZeros( mat  ,  4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  4 ||
          sm(1,0) != -3 || sm(1,1) !=  5 ||
          sm(2,0) !=  0 || sm(2,1) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  4 )\n( -3  5 )\n(  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  0 || mat(0,2) !=  4 || mat(0,3) != 0 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != -3 || mat(1,2) !=  5 || mat(1,3) != 0 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) !=  0 || mat(2,2) != -6 || mat(2,3) != 0 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) !=  0 || mat(3,3) != 0 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  4  0  0 )\n"
                                     "( 0 -3  5  0  0 )\n"
                                     "( 0  0 -6  0  0 )\n"
                                     "( 0  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix copy assignment (aliasing)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );
      sm = blaze::submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  4 ||
          sm(1,0) != -3 || sm(1,1) !=  5 ||
          sm(2,0) !=  0 || sm(2,1) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  4 )\n( -3  5 )\n(  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -3 || tmat_(1,2) !=  5 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -6 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  4  0  7 )\n"
                                     "( 0 -3  5  4 -8 )\n"
                                     "( 0  0 -6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix assignment (mixed type)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {  0, 12 },
                                                      { 11, 13 },
                                                      {  0, 14 } };

      sm = mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 ||
          sm(1,0) != 11 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12 )\n( 11 13 )\n(  0 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 12 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 14 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 12  0  7 )\n"
                                     "( 0 11 13  4 -8 )\n"
                                     "( 0  0 14  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix assignment (aligned/padded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 48UL ), 3UL, 2UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm = mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 ||
          sm(1,0) != 11 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12 )\n( 11 13 )\n(  0 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 12 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 14 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 12  0  7 )\n"
                                     "( 0 11 13  4 -8 )\n"
                                     "( 0  0 14  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 3UL, 2UL );
      mat = 0;
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm = mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 ||
          sm(1,0) != 11 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12 )\n( 11 13 )\n(  0 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 12 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 14 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 12  0  7 )\n"
                                     "( 0 11 13  4 -8 )\n"
                                     "( 0  0 14  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (mixed type)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {  0, 12 },
                                                         { 11, 13 },
                                                         {  0, 14 } };

      sm = mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 ||
          sm(1,0) != 11 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12 )\n( 11 13 )\n(  0 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 12 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 14 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 12  0  7 )\n"
                                     "( 0 11 13  4 -8 )\n"
                                     "( 0  0 14  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (aligned/padded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 3UL, 2UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm = mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 ||
          sm(1,0) != 11 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12 )\n( 11 13 )\n(  0 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 12 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 14 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 12  0  7 )\n"
                                     "( 0 11 13  4 -8 )\n"
                                     "( 0  0 14  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 3UL, 2UL );
      mat = 0;
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm = mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 ||
          sm(1,0) != 11 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12 )\n( 11 13 )\n(  0 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 12 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 14 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 12  0  7 )\n"
                                     "( 0 11 13  4 -8 )\n"
                                     "( 0  0 14  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix assignment";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm = mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 ||
          sm(1,0) != 11 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12 )\n( 11 13 )\n(  0 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 12 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 14 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 12  0  7 )\n"
                                     "( 0 11 13  4 -8 )\n"
                                     "( 0  0 14  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix assignment";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm = mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 ||
          sm(1,0) != 11 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12 )\n( 11 13 )\n(  0 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 12 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 14 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 12  0  7 )\n"
                                     "( 0 11 13  4 -8 )\n"
                                     "( 0  0 14  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testAddAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Submatrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix addition assignment (no aliasing)";

      initialize();

      MT mat( 5UL, 4UL, 0 );
      mat(1,0) = 11;
      mat(2,0) = 12;
      mat(2,2) = 13;

      SMT sm = blaze::submatrix( mat, 1UL, 0UL, 2UL, 3UL );
      sm += blaze::submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  5UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );
      checkRows    ( mat ,  5UL );
      checkColumns ( mat ,  4UL );
      checkNonZeros( mat ,  5UL );

      if( sm(0,0) != 11 || sm(0,1) != -3 || sm(0,2) != 0 ||
          sm(1,0) != 16 || sm(1,1) !=  5 || sm(1,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 11 -3  0 )\n( 16  5  7 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=  0 || mat(0,1) !=  0 || mat(0,2) != 0 || mat(0,3) != 0 ||
          mat(1,0) != 11 || mat(1,1) != -3 || mat(1,2) != 0 || mat(1,3) != 0 ||
          mat(2,0) != 16 || mat(2,1) !=  5 || mat(2,2) != 7 || mat(2,3) != 0 ||
          mat(3,0) !=  0 || mat(3,1) !=  0 || mat(3,2) != 0 || mat(3,3) != 0 ||
          mat(4,0) !=  0 || mat(4,1) !=  0 || mat(4,2) != 0 || mat(4,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 11 -3  0  0 )\n"
                                     "( 16  5  7  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "(  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix addition assignment (aliasing)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );
      sm += blaze::submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) !=  0 ||
          sm(1,0) != 2 || sm(1,1) !=  5 || sm(1,2) != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -2  0 )\n( 2  5 -9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 2 || mat_(2,1) !=  5 || mat_(2,2) != -9 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0 -2  0  0 )\n"
                                     "( 2  5 -9  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix addition assignment (mixed type)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {  0, 11,  0 },
                                                      { 12, 13, 14 } };

      sm += mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix addition assignment (aligned/padded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm += mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 3UL );
      mat = 0;
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm += mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (mixed type)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {  0, 11,  0 },
                                                         { 12, 13, 14 } };

      sm += mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (aligned/padded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm += mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 3UL );
      mat = 0;
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm += mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix addition assignment";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm += mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix addition assignment";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      sm += mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Submatrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix addition assignment (no aliasing)";

      initialize();

      OMT mat( 4UL, 5UL, 0 );
      mat(0,1) = 11;
      mat(0,2) = 12;
      mat(2,2) = 13;

      OSMT sm = blaze::submatrix( mat, 0UL, 1UL, 3UL, 2UL );
      sm += blaze::submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  5UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );
      checkRows    ( mat  ,  4UL );
      checkColumns ( mat  ,  5UL );
      checkNonZeros( mat  ,  5UL );

      if( sm(0,0) != 11 || sm(0,1) != 16 ||
          sm(1,0) != -3 || sm(1,1) !=  5 ||
          sm(2,0) !=  0 || sm(2,1) !=  7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 11 16 )\n( -3  5 )\n(  0  7 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 11 || mat(0,2) != 16 || mat(0,3) != 0 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != -3 || mat(1,2) !=  5 || mat(1,3) != 0 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) !=  0 || mat(2,2) !=  7 || mat(2,3) != 0 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) !=  0 || mat(3,3) != 0 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 11 16  0  0 )\n"
                                     "( 0 -3  5  0  0 )\n"
                                     "( 0  0  7  0  0 )\n"
                                     "( 0  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix addition assignment (aliasing)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );
      sm += blaze::submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  2 ||
          sm(1,0) != -2 || sm(1,1) !=  5 ||
          sm(2,0) !=  0 || sm(2,1) != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  2 )\n( -2  5 )\n(  0 -9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -2 || tmat_(1,2) !=  5 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -9 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  2  0  7 )\n"
                                     "( 0 -2  5  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix addition assignment (mixed type)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {  0, 12 },
                                                      { 11, 13 },
                                                      {  0, 14 } };

      sm += mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix addition assignment (aligned/padded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 48UL ), 3UL, 2UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm += mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 3UL, 2UL );
      mat = 0;
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm += mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (mixed type)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {  0, 12 },
                                                         { 11, 13 },
                                                         {  0, 14 } };

      sm += mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (aligned/padded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 3UL, 2UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm += mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 3UL, 2UL );
      mat = 0;
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm += mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix addition assignment";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm += mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix addition assignment";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      sm += mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testSubAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Submatrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix subtraction assignment (no aliasing)";

      initialize();

      MT mat( 5UL, 4UL, 0 );
      mat(1,0) = 11;
      mat(2,0) = 12;
      mat(2,2) = 13;

      SMT sm = blaze::submatrix( mat, 1UL, 0UL, 2UL, 3UL );
      sm -= blaze::submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  5UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );
      checkRows    ( mat ,  5UL );
      checkColumns ( mat ,  4UL );
      checkNonZeros( mat ,  5UL );

      if( sm(0,0) != 11 || sm(0,1) !=  3 || sm(0,2) !=  0 ||
          sm(1,0) !=  8 || sm(1,1) != -5 || sm(1,2) != 19 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 11  3  0 )\n(  8 -5 19 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=  0 || mat(0,1) !=  0 || mat(0,2) !=  0 || mat(0,3) != 0 ||
          mat(1,0) != 11 || mat(1,1) !=  3 || mat(1,2) !=  0 || mat(1,3) != 0 ||
          mat(2,0) !=  8 || mat(2,1) != -5 || mat(2,2) != 19 || mat(2,3) != 0 ||
          mat(3,0) !=  0 || mat(3,1) !=  0 || mat(3,2) !=  0 || mat(3,3) != 0 ||
          mat(4,0) !=  0 || mat(4,1) !=  0 || mat(4,2) !=  0 || mat(4,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 11  3  0  0 )\n"
                                     "(  8 -5 19  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "(  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix subtraction assignment (aliasing)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );
      sm -= blaze::submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  4 || sm(0,2) != 0 ||
          sm(1,0) != -6 || sm(1,1) != -5 || sm(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  4  0 )\n( -6 -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  4 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -6 || mat_(2,1) != -5 || mat_(2,2) != 3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  4  0  0 )\n"
                                     "( -6 -5  3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (mixed type)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {   0, -11,   0 },
                                                      { -12, -13, -14 } };

      sm -= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      sm -= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 3UL );
      mat = 0;
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      sm -= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (mixed type)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {   0, -11,   0 },
                                                         { -12, -13, -14 } };

      sm -= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 48UL ), 2UL, 3UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      sm -= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 3UL );
      mat = 0;
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      sm -= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 10 13 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix subtraction assignment";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      sm -= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix subtraction assignment";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      sm -= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 ||
          sm(1,0) != 10 || sm(1,1) != 13 || sm(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0 )\n( 12 13 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 12 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 10 || mat_(2,1) != 13 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 10 13 11  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Submatrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix subtraction assignment (no aliasing)";

      initialize();

      OMT mat( 4UL, 5UL, 0 );
      mat(0,1) = 11;
      mat(0,2) = 12;
      mat(2,2) = 13;

      OSMT sm = blaze::submatrix( mat, 0UL, 1UL, 3UL, 2UL );
      sm -= blaze::submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  5UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );
      checkRows    ( mat  ,  4UL );
      checkColumns ( mat  ,  5UL );
      checkNonZeros( mat  ,  5UL );

      if( sm(0,0) != 11 || sm(0,1) !=  8 ||
          sm(1,0) !=  3 || sm(1,1) != -5 ||
          sm(2,0) !=  0 || sm(2,1) != 19 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 11  8 )\n(  3 -5 )\n(  0 19 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 11 || mat(0,2) !=  8 || mat(0,3) != 0 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) !=  3 || mat(1,2) != -5 || mat(1,3) != 0 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) !=  0 || mat(2,2) != 19 || mat(2,3) != 0 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) !=  0 || mat(3,3) != 0 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 11  8  0  0 )\n"
                                     "( 0  3 -5  0  0 )\n"
                                     "( 0  0 19  0  0 )\n"
                                     "( 0  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix subtraction assignment (aliasing)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );
      sm -= blaze::submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) != 0 || sm(0,1) != -6 ||
          sm(1,0) != 4 || sm(1,1) != -5 ||
          sm(2,0) != 0 || sm(2,1) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -6 )\n( 4 -5 )\n( 0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 4 || tmat_(1,2) != -5 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  4 -5  4 -8 )\n"
                                     "( 0  0  3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (mixed type)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {   0, -12 },
                                                      { -11, -13 },
                                                      {   0, -14 } };

      sm -= mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 48UL ), 3UL, 2UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      sm -= mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 3UL, 2UL );
      mat = 0;
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      sm -= mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (mixed type)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {   0, -12 },
                                                         { -11, -13 },
                                                         {   0, -14 } };

      sm -= mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 3UL, 2UL, 16UL, blaze::Deallocate() );
      mat = 0;
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      sm -= mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[7UL] );
      UnalignedUnpadded mat( array.get()+1UL, 3UL, 2UL );
      mat = 0;
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      sm -= mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix subtraction assignment";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      sm -= mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix subtraction assignment";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      sm -= mat;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) !=  0 || sm(0,1) != 10 ||
          sm(1,0) != 12 || sm(1,1) != 13 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 10 )\n( 12 13 )\n(  0 11 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 12 || tmat_(1,2) != 13 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 10  0  7 )\n"
                                     "( 0 12 13  4 -8 )\n"
                                     "( 0  0 11  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testMultAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Submatrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix multiplication assignment (no aliasing)";

      initialize();

      MT mat( 5UL, 4UL, 0 );
      mat(1,0) = 1;
      mat(1,1) = 1;
      mat(2,0) = 1;
      mat(2,1) = 1;

      SMT sm = blaze::submatrix( mat, 1UL, 0UL, 2UL, 2UL );
      sm *= blaze::submatrix( mat_, 2UL, 1UL, 2UL, 2UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );
      checkRows    ( mat ,  5UL );
      checkColumns ( mat ,  4UL );
      checkNonZeros( mat ,  4UL );

      if( sm(0,0) != 4 || sm(0,1) != 2 ||
          sm(1,0) != 4 || sm(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 4 2 )\n( 4 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 2 || mat(1,2) != 0 || mat(1,3) != 0 ||
          mat(2,0) != 4 || mat(2,1) != 2 || mat(2,2) != 0 || mat(2,3) != 0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != 0 || mat(3,3) != 0 ||
          mat(4,0) != 0 || mat(4,1) != 0 || mat(4,2) != 0 || mat(4,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 4  2  0  0 )\n"
                                     "( 4  2  0  0 )\n"
                                     "( 0  0  0  0 )\n"
                                     "( 0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix multiplication assignment (aliasing)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );
      sm *= blaze::submatrix( mat_, 2UL, 1UL, 2UL, 2UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  3UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) != 4 || sm(0,1) != 5 ||
          sm(1,0) != 0 || sm(1,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 4  5 )\n( 0  6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 4 || mat_(1,1) !=  5 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) !=  6 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 4  5  0  0 )\n"
                                     "( 0  6 -3  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (mixed type)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { -11, -12 },
                                                      {  13,  14 } };

      sm *= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 13 || sm(0,1) != 14 ||
          sm(1,0) != 22 || sm(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 13 14 )\n( 22 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 13 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 22 || mat_(2,1) != 24 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 13 14  0  0 )\n"
                                     "( 22 24 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 2UL, 2UL, 16UL, blaze::Deallocate() );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      sm *= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 13 || sm(0,1) != 14 ||
          sm(1,0) != 22 || sm(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 13 14 )\n( 22 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 13 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 22 || mat_(2,1) != 24 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 13 14  0  0 )\n"
                                     "( 22 24 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[5UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 2UL );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      sm *= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 13 || sm(0,1) != 14 ||
          sm(1,0) != 22 || sm(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 13 14 )\n( 22 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 13 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 22 || mat_(2,1) != 24 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 13 14  0  0 )\n"
                                     "( 22 24 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (mixed type)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { -11, -12 },
                                                         {  13,  14 } };

      sm *= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 13 || sm(0,1) != 14 ||
          sm(1,0) != 22 || sm(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 13 14 )\n( 22 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 13 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 22 || mat_(2,1) != 24 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 13 14  0  0 )\n"
                                     "( 22 24 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 2UL, 2UL, 16UL, blaze::Deallocate() );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      sm *= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 13 || sm(0,1) != 14 ||
          sm(1,0) != 22 || sm(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 13 14 )\n( 22 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 13 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 22 || mat_(2,1) != 24 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 13 14  0  0 )\n"
                                     "( 22 24 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[5UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 2UL );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      sm *= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 13 || sm(0,1) != 14 ||
          sm(1,0) != 22 || sm(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 13 14 )\n( 22 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 13 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 22 || mat_(2,1) != 24 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 13 14  0  0 )\n"
                                     "( 22 24 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix multiplication assignment";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      sm *= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 13 || sm(0,1) != 14 ||
          sm(1,0) != 22 || sm(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 13 14 )\n( 22 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 13 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 22 || mat_(2,1) != 24 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 13 14  0  0 )\n"
                                     "( 22 24 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix multiplication assignment";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      sm *= mat;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 13 || sm(0,1) != 14 ||
          sm(1,0) != 22 || sm(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 13 14 )\n( 22 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 13 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 22 || mat_(2,1) != 24 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 13 14  0  0 )\n"
                                     "( 22 24 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Submatrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix multiplication assignment (no aliasing)";

      initialize();

      OMT mat( 4UL, 5UL, 0 );
      mat(0,1) = 1;
      mat(0,2) = 1;
      mat(1,1) = 1;
      mat(1,2) = 1;

      OSMT sm = blaze::submatrix( mat, 0UL, 1UL, 2UL, 2UL );
      sm *= blaze::submatrix( tmat_, 1UL, 2UL, 2UL, 2UL );

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );
      checkRows    ( mat  ,  4UL );
      checkColumns ( mat  ,  5UL );
      checkNonZeros( mat  ,  4UL );

      if( sm(0,0) != -3 || sm(0,1) != 9 ||
          sm(1,0) != -3 || sm(1,1) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -3 -3 )\n(  9  9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != -3 || mat(0,2) != 9 || mat(0,3) != 0 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != -3 || mat(1,2) != 9 || mat(1,3) != 0 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) !=  0 || mat(2,2) != 0 || mat(2,3) != 0 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) != 0 || mat(3,3) != 0 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 -3  9  0  0 )\n"
                                     "( 0 -3  9  0  0 )\n"
                                     "( 0  0  0  0  0 )\n"
                                     "( 0  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix multiplication assignment (aliasing)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );
      sm *= blaze::submatrix( tmat_, 1UL, 2UL, 2UL, 2UL );

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  3UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) != 6 || sm(0,1) != -10 ||
          sm(1,0) != 0 || sm(1,1) !=   4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 6 -10 )\n( 0   4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 6 || tmat_(0,2) != -10 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 0 || tmat_(1,2) !=   4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=   0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  6 -10  0  7 )\n"
                                     "( 0  0   4  4 -8 )\n"
                                     "( 0  0  -3  5  9 )\n"
                                     "( 0  0   0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (mixed type)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {  11,  12 },
                                                      { -13, -14 } };

      sm *= mat;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 26 || sm(0,1) != 28 ||
          sm(1,0) != 11 || sm(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 26 28 )\n( 11 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 26 || tmat_(0,2) != 28 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 12 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 26 28  0  7 )\n"
                                     "( 0 11 12  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,rowMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 2UL, 2UL, 16UL, blaze::Deallocate() );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      sm *= mat;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 26 || sm(0,1) != 28 ||
          sm(1,0) != 11 || sm(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 26 28 )\n( 11 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 26 || tmat_(0,2) != 28 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 12 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 26 28  0  7 )\n"
                                     "( 0 11 12  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[5UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 2UL );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      sm *= mat;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 26 || sm(0,1) != 28 ||
          sm(1,0) != 11 || sm(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 26 28 )\n( 11 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 26 || tmat_(0,2) != 28 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 12 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 26 28  0  7 )\n"
                                     "( 0 11 12  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (mixed type)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {  11,  12 },
                                                         { -13, -14 } };

      sm *= mat;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 26 || sm(0,1) != 28 ||
          sm(1,0) != 11 || sm(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 26 28 )\n( 11 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 26 || tmat_(0,2) != 28 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 12 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 26 28  0  7 )\n"
                                     "( 0 11 12  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      typedef blaze::CustomMatrix<int,aligned,padded,columnMajor>  AlignedPadded;
      AlignedPadded mat( blaze::allocate<int>( 32UL ), 2UL, 2UL, 16UL, blaze::Deallocate() );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      sm *= mat;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 26 || sm(0,1) != 28 ||
          sm(1,0) != 11 || sm(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 26 28 )\n( 11 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 26 || tmat_(0,2) != 28 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 12 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 26 28  0  7 )\n"
                                     "( 0 11 12  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (unaligned/unpadded))";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      typedef blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<int[]> array( new int[5UL] );
      UnalignedUnpadded mat( array.get()+1UL, 2UL, 2UL );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      sm *= mat;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 26 || sm(0,1) != 28 ||
          sm(1,0) != 11 || sm(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 26 28 )\n( 11 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 26 || tmat_(0,2) != 28 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 12 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 26 28  0  7 )\n"
                                     "( 0 11 12  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix multiplication assignment";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      sm *= mat;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 26 || sm(0,1) != 28 ||
          sm(1,0) != 11 || sm(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 26 28 )\n( 11 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 26 || tmat_(0,2) != 28 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 12 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 26 28  0  7 )\n"
                                     "( 0 11 12  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix multiplication assignment";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      sm *= mat;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 26 || sm(0,1) != 28 ||
          sm(1,0) != 11 || sm(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 26 28 )\n( 11 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 26 || tmat_(0,2) != 28 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 11 || tmat_(1,2) != 12 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 26 28  0  7 )\n"
                                     "( 0 11 12  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all Submatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s) (2x3)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 2UL, 3UL );

      sm *= 3;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 || sm(0,2) != -9 ||
          sm(1,0) !=  0 || sm(1,1) != 12 || sm(1,2) != 15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0 -9 )\n(  0 12 15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -6 || mat_(2,1) !=  0 || mat_(2,2) != -9 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 12 || mat_(3,2) != 15 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -6  0 -9  0 )\n"
                                     "(  0 12 15 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-scaling (M*=s) (3x2)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 3UL, 2UL );

      sm *= 3;

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=   0 ||
          sm(1,0) !=  0 || sm(1,1) !=  12 ||
          sm(2,0) != 21 || sm(2,1) != -24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6   0 )\n(  0  12 )\n( 21 -24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=   0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=   1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -6 || mat_(2,1) !=   0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  12 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 21 || mat_(4,1) != -24 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0   0  0  0 )\n"
                                     "(  0   1  0  0 )\n"
                                     "( -6   0 -3  0 )\n"
                                     "(  0  12  5 -6 )\n"
                                     "( 21 -24  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s) (2x3)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 2UL, 3UL );

      sm = sm * 3;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 || sm(0,2) != -9 ||
          sm(1,0) !=  0 || sm(1,1) != 12 || sm(1,2) != 15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0 -9 )\n(  0 12 15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -6 || mat_(2,1) !=  0 || mat_(2,2) != -9 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 12 || mat_(3,2) != 15 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -6  0 -9  0 )\n"
                                     "(  0 12 15 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-scaling (M=M*s) (3x2)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 3UL, 2UL );

      sm = sm * 3;

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=   0 ||
          sm(1,0) !=  0 || sm(1,1) !=  12 ||
          sm(2,0) != 21 || sm(2,1) != -24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6   0 )\n(  0  12 )\n( 21 -24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=   0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=   1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -6 || mat_(2,1) !=   0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  12 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 21 || mat_(4,1) != -24 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0   0  0  0 )\n"
                                     "(  0   1  0  0 )\n"
                                     "( -6   0 -3  0 )\n"
                                     "(  0  12  5 -6 )\n"
                                     "( 21 -24  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M) (2x3)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 2UL, 3UL );

      sm = 3 * sm;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 || sm(0,2) != -9 ||
          sm(1,0) !=  0 || sm(1,1) != 12 || sm(1,2) != 15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0 -9 )\n(  0 12 15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -6 || mat_(2,1) !=  0 || mat_(2,2) != -9 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 12 || mat_(3,2) != 15 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -6  0 -9  0 )\n"
                                     "(  0 12 15 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-scaling (M=s*M) (3x2)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 3UL, 2UL );

      sm = 3 * sm;

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=   0 ||
          sm(1,0) !=  0 || sm(1,1) !=  12 ||
          sm(2,0) != 21 || sm(2,1) != -24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6   0 )\n(  0  12 )\n( 21 -24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=   0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=   1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -6 || mat_(2,1) !=   0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  12 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 21 || mat_(4,1) != -24 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0   0  0  0 )\n"
                                     "(  0   1  0  0 )\n"
                                     "( -6   0 -3  0 )\n"
                                     "(  0  12  5 -6 )\n"
                                     "( 21 -24  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s) (2x3)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 2UL, 3UL );

      sm /= 0.5;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -4 || sm(0,1) != 0 || sm(0,2) != -6 ||
          sm(1,0) !=  0 || sm(1,1) != 8 || sm(1,2) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4  0 -6 )\n(  0  8 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=  0 || mat_(2,2) != -6 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  8 || mat_(3,2) != 10 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  0 -6  0 )\n"
                                     "(  0  8 10 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-scaling (M/=s) (3x2)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 3UL, 2UL );

      sm /= 0.5;

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -4 || sm(0,1) !=   0 ||
          sm(1,0) !=  0 || sm(1,1) !=   8 ||
          sm(2,0) != 14 || sm(2,1) != -16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4   0 )\n(  0   8 )\n( 14 -16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=   0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=   1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=   0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=   8 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 14 || mat_(4,1) != -16 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0   0  0  0 )\n"
                                     "(  0   1  0  0 )\n"
                                     "( -4   0 -3  0 )\n"
                                     "(  0   8  5 -6 )\n"
                                     "( 14 -16  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s) (2x3)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 2UL, 3UL );

      sm = sm / 0.5;

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -4 || sm(0,1) != 0 || sm(0,2) != -6 ||
          sm(1,0) !=  0 || sm(1,1) != 8 || sm(1,2) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4  0 -6 )\n(  0  8 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=  0 || mat_(2,2) != -6 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  8 || mat_(3,2) != 10 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -4  0 -6  0 )\n"
                                     "(  0  8 10 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-scaling (M=M/s) (3x2)";

      initialize();

      SMT sm = blaze::submatrix( mat_, 2UL, 0UL, 3UL, 2UL );

      sm = sm / 0.5;

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  2UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != -4 || sm(0,1) !=   0 ||
          sm(1,0) !=  0 || sm(1,1) !=   8 ||
          sm(2,0) != 14 || sm(2,1) != -16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4   0 )\n(  0   8 )\n( 14 -16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=   0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=   1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -4 || mat_(2,1) !=   0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=   8 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) != 14 || mat_(4,1) != -16 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0   0  0  0 )\n"
                                     "(  0   1  0  0 )\n"
                                     "( -4   0 -3  0 )\n"
                                     "(  0   8  5 -6 )\n"
                                     "( 14 -16  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major Submatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major Submatrix::scale()";

      initialize();

      // Initialization check
      SMT sm = blaze::submatrix( mat_, 2UL, 1UL, 2UL, 2UL );

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 3UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 2UL );

      if( sm(0,0) != 0 || sm(0,1) != -3 ||
          sm(1,0) != 4 || sm(1,1) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -3 )\n( 4  5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      sm.scale( 2 );

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 3UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 2UL );

      if( sm(0,0) != 0 || sm(0,1) != -6 ||
          sm(1,0) != 8 || sm(1,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Integral scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -6 )\n( 8 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      sm.scale( 0.5 );

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 3UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 2UL );

      if( sm(0,0) != 0 || sm(0,1) != -3 ||
          sm(1,0) != 4 || sm(1,1) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Floating point scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -3 )\n( 4  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s) (3x2)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 3UL, 2UL );

      sm *= 3;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 ||
          sm(1,0) !=  0 || sm(1,1) != 12 ||
          sm(2,0) != -9 || sm(2,1) != 15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0 )\n(  0 12 )\n( -9 15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) != 12 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -9 || tmat_(2,3) != 15 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0 12 -8 )\n"
                                     "( 0  0 -9 15  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-scaling (M*=s) (2x3)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 2UL, 3UL );

      sm *= 3;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 || sm(0,2) !=  21 ||
          sm(1,0) !=  0 || sm(1,1) != 12 || sm(1,2) != -24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0  21 )\n(  0 12 -24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  21 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) != 12 || tmat_(1,4) != -24 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=   9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  21 )\n"
                                     "( 0  1  0 12 -24 )\n"
                                     "( 0  0 -3  5   9 )\n"
                                     "( 0  0  0 -6  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s) (3x2)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 3UL, 2UL );

      sm = sm * 3;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 ||
          sm(1,0) !=  0 || sm(1,1) != 12 ||
          sm(2,0) != -9 || sm(2,1) != 15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0 )\n(  0 12 )\n( -9 15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) != 12 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -9 || tmat_(2,3) != 15 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0 12 -8 )\n"
                                     "( 0  0 -9 15  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-scaling (M=M*s) (2x3)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 2UL, 3UL );

      sm = sm * 3;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 || sm(0,2) !=  21 ||
          sm(1,0) !=  0 || sm(1,1) != 12 || sm(1,2) != -24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0  21 )\n(  0 12 -24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  21 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) != 12 || tmat_(1,4) != -24 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=   9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  21 )\n"
                                     "( 0  1  0 12 -24 )\n"
                                     "( 0  0 -3  5   9 )\n"
                                     "( 0  0  0 -6  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M) (3x2)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 3UL, 2UL );

      sm = 3 * sm;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 ||
          sm(1,0) !=  0 || sm(1,1) != 12 ||
          sm(2,0) != -9 || sm(2,1) != 15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0 )\n(  0 12 )\n( -9 15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) != 12 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -9 || tmat_(2,3) != 15 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0 12 -8 )\n"
                                     "( 0  0 -9 15  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-scaling (M=s*M) (2x3)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 2UL, 3UL );

      sm = 3 * sm;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -6 || sm(0,1) !=  0 || sm(0,2) !=  21 ||
          sm(1,0) !=  0 || sm(1,1) != 12 || sm(1,2) != -24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -6  0  21 )\n(  0 12 -24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  21 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) != 12 || tmat_(1,4) != -24 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=   9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  21 )\n"
                                     "( 0  1  0 12 -24 )\n"
                                     "( 0  0 -3  5   9 )\n"
                                     "( 0  0  0 -6  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s) (3x2)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 3UL, 2UL );

      sm /= 0.5;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -4 || sm(0,1) !=  0 ||
          sm(1,0) !=  0 || sm(1,1) !=  8 ||
          sm(2,0) != -6 || sm(2,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4  0 )\n(  0  8 )\n( -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -6 || tmat_(2,3) != 10 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  8 -8 )\n"
                                     "( 0  0 -6 10  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-scaling (M/=s) (2x3)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 2UL, 3UL );

      sm /= 0.5;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -4 || sm(0,1) != 0 || sm(0,2) !=  14 ||
          sm(1,0) !=  0 || sm(1,1) != 8 || sm(1,2) != -16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4  0  14 )\n(  0  8 -16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  14 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -16 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=   9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  14 )\n"
                                     "( 0  1  0  8 -16 )\n"
                                     "( 0  0 -3  5   9 )\n"
                                     "( 0  0  0 -6  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s) (3x2)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 3UL, 2UL );

      sm = sm / 0.5;

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -4 || sm(0,1) !=  0 ||
          sm(1,0) !=  0 || sm(1,1) !=  8 ||
          sm(2,0) != -6 || sm(2,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4  0 )\n(  0  8 )\n( -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -6 || tmat_(2,3) != 10 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  8 -8 )\n"
                                     "( 0  0 -6 10  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-scaling (M=M/s) (2x3)";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 2UL, 2UL, 3UL );

      sm = sm / 0.5;

      checkRows    ( sm   ,  2UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) != -4 || sm(0,1) != 0 || sm(0,2) !=  14 ||
          sm(1,0) !=  0 || sm(1,1) != 8 || sm(1,2) != -16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4  0  14 )\n(  0  8 -16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  14 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -16 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=   9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  14 )\n"
                                     "( 0  1  0  8 -16 )\n"
                                     "( 0  0 -3  5   9 )\n"
                                     "( 0  0  0 -6  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Submatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major Submatrix::scale()";

      initialize();

      // Initialization check
      OSMT sm = blaze::submatrix( tmat_, 1UL, 2UL, 2UL, 2UL );

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 3UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 2UL );

      if( sm(0,0) !=  0 || sm(0,1) != 4 ||
          sm(1,0) != -3 || sm(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 4 )\n( -3 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      sm.scale( 2 );

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 3UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 2UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  8 ||
          sm(1,0) != -6 || sm(1,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Integral scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  8 )\n( -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      sm.scale( 0.5 );

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 3UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 2UL );

      if( sm(0,0) !=  0 || sm(0,1) != 4 ||
          sm(1,0) != -3 || sm(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Floating point scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 4 )\n( -3 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the Submatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void DenseUnalignedTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major Submatrix::operator()";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 2UL );

      // Assignment to the element (1,0)
      {
         sm(1,0) = 9;

         checkRows    ( sm  ,  3UL );
         checkColumns ( sm  ,  2UL );
         checkNonZeros( sm  ,  5UL );
         checkNonZeros( sm  ,  0UL, 1UL );
         checkNonZeros( sm  ,  1UL, 2UL );
         checkNonZeros( sm  ,  2UL, 2UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 11UL );

         if( sm(0,0) != 1 || sm(0,1) !=  0 ||
             sm(1,0) != 9 || sm(1,1) != -3 ||
             sm(2,0) != 4 || sm(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 1  0 )\n( 9 -3 )\n( 4  5 )\n";
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
      }

      // Assignment to the element (2,0)
      {
         sm(2,0) = 0;

         checkRows    ( sm  ,  3UL );
         checkColumns ( sm  ,  2UL );
         checkNonZeros( sm  ,  4UL );
         checkNonZeros( sm  ,  0UL, 1UL );
         checkNonZeros( sm  ,  1UL, 2UL );
         checkNonZeros( sm  ,  2UL, 1UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 10UL );

         if( sm(0,0) != 1 || sm(0,1) !=  0 ||
             sm(1,0) != 9 || sm(1,1) != -3 ||
             sm(2,0) != 0 || sm(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 1  0 )\n( 9 -3 )\n( 0  5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  9 -3  0 )\n"
                                        "(  0  0  5 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Assignment to the element (1,1)
      {
         sm(1,1) = 11;

         checkRows    ( sm  ,  3UL );
         checkColumns ( sm  ,  2UL );
         checkNonZeros( sm  ,  4UL );
         checkNonZeros( sm  ,  0UL, 1UL );
         checkNonZeros( sm  ,  1UL, 2UL );
         checkNonZeros( sm  ,  2UL, 1UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 10UL );

         if( sm(0,0) != 1 || sm(0,1) !=  0 ||
             sm(1,0) != 9 || sm(1,1) != 11 ||
             sm(2,0) != 0 || sm(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 1  0 )\n( 9 11 )\n( 0  5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  9 11  0 )\n"
                                        "(  0  0  5 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Addition assignment to the element (0,0)
      {
         sm(0,0) += 3;

         checkRows    ( sm  ,  3UL );
         checkColumns ( sm  ,  2UL );
         checkNonZeros( sm  ,  4UL );
         checkNonZeros( sm  ,  0UL, 1UL );
         checkNonZeros( sm  ,  1UL, 2UL );
         checkNonZeros( sm  ,  2UL, 1UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 10UL );

         if( sm(0,0) != 4 || sm(0,1) !=  0 ||
             sm(1,0) != 9 || sm(1,1) != 11 ||
             sm(2,0) != 0 || sm(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 4  0 )\n( 9 11 )\n( 0  5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  4 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  4  0  0 )\n"
                                        "( -2  9 11  0 )\n"
                                        "(  0  0  5 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Subtraction assignment to the element (0,1)
      {
         sm(0,1) -= 6;

         checkRows    ( sm  ,  3UL );
         checkColumns ( sm  ,  2UL );
         checkNonZeros( sm  ,  5UL );
         checkNonZeros( sm  ,  0UL, 2UL );
         checkNonZeros( sm  ,  1UL, 2UL );
         checkNonZeros( sm  ,  2UL, 1UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 11UL );

         if( sm(0,0) != 4 || sm(0,1) != -6 ||
             sm(1,0) != 9 || sm(1,1) != 11 ||
             sm(2,0) != 0 || sm(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 4 -6 )\n( 9 11 )\n( 0  5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  4 || mat_(1,2) != -6 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  4 -6  0 )\n"
                                        "( -2  9 11  0 )\n"
                                        "(  0  0  5 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Multiplication assignment to the element (1,1)
      {
         sm(1,1) *= 2;

         checkRows    ( sm  ,  3UL );
         checkColumns ( sm  ,  2UL );
         checkNonZeros( sm  ,  5UL );
         checkNonZeros( sm  ,  0UL, 2UL );
         checkNonZeros( sm  ,  1UL, 2UL );
         checkNonZeros( sm  ,  2UL, 1UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 11UL );

         if( sm(0,0) != 4 || sm(0,1) != -6 ||
             sm(1,0) != 9 || sm(1,1) != 22 ||
             sm(2,0) != 0 || sm(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 4 -6 )\n( 9 22 )\n( 0  5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  4 || mat_(1,2) != -6 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != 22 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  4 -6  0 )\n"
                                        "( -2  9 22  0 )\n"
                                        "(  0  0  5 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Division assignment to the element (1,1)
      {
         sm(1,1) /= 2;

         checkRows    ( sm  ,  3UL );
         checkColumns ( sm  ,  2UL );
         checkNonZeros( sm  ,  5UL );
         checkNonZeros( sm  ,  0UL, 2UL );
         checkNonZeros( sm  ,  1UL, 2UL );
         checkNonZeros( sm  ,  2UL, 1UL );
         checkRows    ( mat_,  5UL );
         checkColumns ( mat_,  4UL );
         checkNonZeros( mat_, 11UL );

         if( sm(0,0) != 4 || sm(0,1) != -6 ||
             sm(1,0) != 9 || sm(1,1) != 11 ||
             sm(2,0) != 0 || sm(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 4 -6 )\n( 9 11 )\n( 0  5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  4 || mat_(1,2) != -6 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  9 || mat_(2,2) != 11 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  4 -6  0 )\n"
                                        "( -2  9 11  0 )\n"
                                        "(  0  0  5 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major Submatrix::operator()";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 1UL, 1UL, 2UL, 3UL );

      // Assignment to the element (0,1)
      {
         sm(0,1) = 9;

         checkRows    ( sm   ,  2UL );
         checkColumns ( sm   ,  3UL );
         checkNonZeros( sm   ,  5UL );
         checkNonZeros( sm   ,  0UL, 1UL );
         checkNonZeros( sm   ,  1UL, 2UL );
         checkNonZeros( sm   ,  2UL, 2UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 11UL );

         if( sm(0,0) != 1 || sm(0,1) !=  9 || sm(0,2) != 4 ||
             sm(1,0) != 0 || sm(1,1) != -3 || sm(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 1  9 4 )\n( 0 -3 5 )\n";
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
      }

      // Assignment to the element (0,2)
      {
         sm(0,2) = 0;

         checkRows    ( sm   ,  2UL );
         checkColumns ( sm   ,  3UL );
         checkNonZeros( sm   ,  4UL );
         checkNonZeros( sm   ,  0UL, 1UL );
         checkNonZeros( sm   ,  1UL, 2UL );
         checkNonZeros( sm   ,  2UL, 1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 10UL );

         if( sm(0,0) != 1 || sm(0,1) !=  9 || sm(0,2) != 0 ||
             sm(1,0) != 0 || sm(1,1) != -3 || sm(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 1  9 0 )\n( 0 -3 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  1  9  0 -8 )\n"
                                        "( 0  0 -3  5  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Assignment to the element (1,1)
      {
         sm(1,1) = 11;

         checkRows    ( sm   ,  2UL );
         checkColumns ( sm   ,  3UL );
         checkNonZeros( sm   ,  4UL );
         checkNonZeros( sm   ,  0UL, 1UL );
         checkNonZeros( sm   ,  1UL, 2UL );
         checkNonZeros( sm   ,  2UL, 1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 10UL );

         if( sm(0,0) != 1 || sm(0,1) !=  9 || sm(0,2) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 11 || sm(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 1 11 0 )\n( 0 -3 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  1  9  0 -8 )\n"
                                        "( 0  0 11  5  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Addition assignment to the element (0,0)
      {
         sm(0,0) += 3;

         checkRows    ( sm   ,  2UL );
         checkColumns ( sm   ,  3UL );
         checkNonZeros( sm   ,  4UL );
         checkNonZeros( sm   ,  0UL, 1UL );
         checkNonZeros( sm   ,  1UL, 2UL );
         checkNonZeros( sm   ,  2UL, 1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 10UL );

         if( sm(0,0) != 4 || sm(0,1) !=  9 || sm(0,2) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 11 || sm(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 4 11 0 )\n( 0 -3 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 4 || tmat_(1,2) !=  9 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  4  9  0 -8 )\n"
                                        "( 0  0 11  5  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Subtraction assignment to the element (1,0)
      {
         sm(1,0) -= 6;

         checkRows    ( sm   ,  2UL );
         checkColumns ( sm   ,  3UL );
         checkNonZeros( sm   ,  5UL );
         checkNonZeros( sm   ,  0UL, 2UL );
         checkNonZeros( sm   ,  1UL, 2UL );
         checkNonZeros( sm   ,  2UL, 1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 11UL );

         if( sm(0,0) !=  4 || sm(0,1) !=  9 || sm(0,2) != 0 ||
             sm(1,0) != -6 || sm(1,1) != 11 || sm(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  4 11 0 )\n( -6 -3 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  4 || tmat_(1,2) !=  9 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != -6 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  4  9  0 -8 )\n"
                                        "( 0 -6 11  5  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Multiplication assignment to the element (1,1)
      {
         sm(1,1) *= 2;

         checkRows    ( sm   ,  2UL );
         checkColumns ( sm   ,  3UL );
         checkNonZeros( sm   ,  5UL );
         checkNonZeros( sm   ,  0UL, 2UL );
         checkNonZeros( sm   ,  1UL, 2UL );
         checkNonZeros( sm   ,  2UL, 1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 11UL );

         if( sm(0,0) !=  4 || sm(0,1) !=  9 || sm(0,2) != 0 ||
             sm(1,0) != -6 || sm(1,1) != 22 || sm(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  4 22 0 )\n( -6 -3 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  4 || tmat_(1,2) !=  9 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != -6 || tmat_(2,2) != 22 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  4  9  0 -8 )\n"
                                        "( 0 -6 22  5  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Division assignment to the element (1,1)
      {
         sm(1,1) /= 2;

         checkRows    ( sm   ,  2UL );
         checkColumns ( sm   ,  3UL );
         checkNonZeros( sm   ,  5UL );
         checkNonZeros( sm   ,  0UL, 2UL );
         checkNonZeros( sm   ,  1UL, 2UL );
         checkNonZeros( sm   ,  2UL, 1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 11UL );

         if( sm(0,0) !=  4 || sm(0,1) !=  9 || sm(0,2) != 0 ||
             sm(1,0) != -6 || sm(1,1) != 11 || sm(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  4 11 0 )\n( -6 -3 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  4 || tmat_(1,2) !=  9 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != -6 || tmat_(2,2) != 11 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  4  9  0 -8 )\n"
                                        "( 0 -6 11  5  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testIterator()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 3UL );

      // Testing the Iterator default constructor
      {
         test_ = "Row-major Iterator default constructor";

         SMT::Iterator it = SMT::Iterator();

         if( it != SMT::Iterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing the ConstIterator default constructor
      {
         test_ = "Row-major ConstIterator default constructor";

         SMT::ConstIterator it = SMT::ConstIterator();

         if( it != SMT::ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Row-major Iterator/ConstIterator conversion";

         SMT::ConstIterator it( begin( sm, 1UL ) );

         if( it == end( sm, 1UL ) || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator
      {
         test_ = "Row-major Iterator subtraction";

         const size_t number( end( sm, 0UL ) - begin( sm, 0UL ) );

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

         const size_t number( cend( sm, 1UL ) - cbegin( sm, 1UL ) );

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

         SMT::ConstIterator it ( cbegin( sm, 2UL ) );
         SMT::ConstIterator end( cend( sm, 2UL ) );

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

         for( SMT::Iterator it=begin( sm, 2UL ); it!=end( sm, 2UL ); ++it ) {
            *it = value++;
         }

         if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
             sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -3 ||
             sm(2,0) !=  7 || sm(2,1) != 8 || sm(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  7 || mat_(3,1) !=  8 || mat_(3,2) !=  9 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  7  8  9 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         int value = 4;

         for( SMT::Iterator it=begin( sm, 1UL ); it!=end( sm, 1UL ); ++it ) {
            *it += value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 0 ||
             sm(1,0) != 2 || sm(1,1) != 5 || sm(1,2) != 3 ||
             sm(2,0) != 7 || sm(2,1) != 8 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 2 5 3 )\n( 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
             mat_(1,0) != 0 || mat_(1,1) !=  1 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
             mat_(2,0) != 2 || mat_(2,1) !=  5 || mat_(2,2) != 3 || mat_(2,3) !=  0 ||
             mat_(3,0) != 7 || mat_(3,1) !=  8 || mat_(3,2) != 9 || mat_(3,3) != -6 ||
             mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "(  2  5  3  0 )\n"
                                        "(  7  8  9 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         int value = 4;

         for( SMT::Iterator it=begin( sm, 1UL ); it!=end( sm, 1UL ); ++it ) {
            *it -= value++;
         }

         if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
             sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -3 ||
             sm(2,0) !=  7 || sm(2,1) != 8 || sm(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  7 || mat_(3,1) !=  8 || mat_(3,2) !=  9 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  7  8  9 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         int value = 2;

         for( SMT::Iterator it=begin( sm, 1UL ); it!=end( sm, 1UL ); ++it ) {
            *it *= value++;
         }

         if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=   0 ||
             sm(1,0) != -4 || sm(1,1) != 0 || sm(1,2) != -12 ||
             sm(2,0) !=  7 || sm(2,1) != 8 || sm(2,2) !=   9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  1   0 )\n( -4  0 -12 )\n(  7  8   9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=   0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -4 || mat_(2,1) !=  0 || mat_(2,2) != -12 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  7 || mat_(3,1) !=  8 || mat_(3,2) !=   9 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0   0  0 )\n"
                                        "(  0  1   0  0 )\n"
                                        "( -4  0 -12  0 )\n"
                                        "(  7  8   9 -6 )\n"
                                        "(  7 -8   9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         for( SMT::Iterator it=begin( sm, 1UL ); it!=end( sm, 1UL ); ++it ) {
            *it /= 2;
         }

         if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
             sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -6 ||
             sm(2,0) !=  7 || sm(2,1) != 8 || sm(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -6 )\n(  7  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -6 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  7 || mat_(3,1) !=  8 || mat_(3,2) !=  9 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -6  0 )\n"
                                        "(  7  8  9 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 3UL );

      // Testing the Iterator default constructor
      {
         test_ = "Column-major Iterator default constructor";

         OSMT::Iterator it = OSMT::Iterator();

         if( it != OSMT::Iterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing the ConstIterator default constructor
      {
         test_ = "Column-major ConstIterator default constructor";

         OSMT::ConstIterator it = OSMT::ConstIterator();

         if( it != OSMT::ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Column-major Iterator/ConstIterator conversion";

         OSMT::ConstIterator it( begin( sm, 1UL ) );

         if( it == end( sm, 1UL ) || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator
      {
         test_ = "Column-major Iterator subtraction";

         const size_t number( end( sm, 0UL ) - begin( sm, 0UL ) );

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

         const size_t number( cend( sm, 1UL ) - cbegin( sm, 1UL ) );

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

         OSMT::ConstIterator it ( cbegin( sm, 2UL ) );
         OSMT::ConstIterator end( cend( sm, 2UL ) );

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

         for( OSMT::Iterator it=begin( sm, 2UL ); it!=end( sm, 2UL ); ++it ) {
            *it = value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 7 ||
             sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -2  7 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  7 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  7  7 )\n"
                                        "( 0  1  0  8 -8 )\n"
                                        "( 0  0 -3  9  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         int value = 4;

         for( OSMT::Iterator it=begin( sm, 1UL ); it!=end( sm, 1UL ); ++it ) {
            *it += value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != 2 || sm(0,2) != 7 ||
             sm(1,0) != 1 || sm(1,1) != 5 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != 3 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 2 7 )\n( 1 5 8 )\n( 0 3 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 2 || tmat_(0,3) !=  7 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) != 5 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 3 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0  2  7  7 )\n"
                                        "( 0  1  5  8 -8 )\n"
                                        "( 0  0  3  9  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         int value = 4;

         for( OSMT::Iterator it=begin( sm, 1UL ); it!=end( sm, 1UL ); ++it ) {
            *it -= value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 7 ||
             sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -2  7 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  7 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  7  7 )\n"
                                        "( 0  1  0  8 -8 )\n"
                                        "( 0  0 -3  9  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         int value = 2;

         for( OSMT::Iterator it=begin( sm, 1UL ); it!=end( sm, 1UL ); ++it ) {
            *it *= value++;
         }

         if( sm(0,0) != 0 || sm(0,1) !=  -4 || sm(0,2) != 7 ||
             sm(1,0) != 1 || sm(1,1) !=   0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != -12 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -2  7 )\n( 1  0  8 )\n( 0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  -4 || tmat_(0,3) !=  7 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=   0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -12 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=   0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0  -4  7  7 )\n"
                                        "( 0  1   0  8 -8 )\n"
                                        "( 0  0 -12  9  9 )\n"
                                        "( 0  0   0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         for( OSMT::Iterator it=begin( sm, 1UL ); it!=end( sm, 1UL ); ++it ) {
            *it /= 2;
         }

         if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 7 ||
             sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != -6 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -2  7 )\n( 1  0  8 )\n( 0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  7 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -6 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  7  7 )\n"
                                        "( 0  1  0  8 -8 )\n"
                                        "( 0  0 -6  9  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testNonZeros()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major Submatrix::nonZeros()";

      initialize();

      // Initialization check
      SMT sm = blaze::submatrix( mat_, 1UL, 1UL, 2UL, 3UL );

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 3UL );
      checkNonZeros( sm, 2UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 1UL );

      if( sm(0,0) != 1 || sm(0,1) !=  0 || sm(0,2) != 0 ||
          sm(1,0) != 0 || sm(1,1) != -3 || sm(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1  0 0 )\n( 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse submatrix
      sm(1,1) = 0;

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 3UL );
      checkNonZeros( sm, 1UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 0UL );

      if( sm(0,0) != 1 || sm(0,1) != 0 || sm(0,2) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      mat_(2,3) = 5;

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 3UL );
      checkNonZeros( sm, 2UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 1UL );

      if( sm(0,0) != 1 || sm(0,1) != 0 || sm(0,2) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major Submatrix::nonZeros()";

      initialize();

      // Initialization check
      OSMT sm = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );

      checkRows    ( sm, 3UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 2UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 1UL );

      if( sm(0,0) != 1 || sm(0,1) !=  0 ||
          sm(1,0) != 0 || sm(1,1) != -3 ||
          sm(2,0) != 0 || sm(2,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1  0 )\n( 0 -3 )\n( 0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse submatrix
      sm(1,1) = 0;

      checkRows    ( sm, 3UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 1UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 0UL );

      if( sm(0,0) != 1 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ||
          sm(2,0) != 0 || sm(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      tmat_(3,2) = 5;

      checkRows    ( sm, 3UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 2UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 1UL );

      if( sm(0,0) != 1 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ||
          sm(2,0) != 0 || sm(2,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 0 )\n( 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testReset()
{
   //=====================================================================================
   // Row-major single element reset
   //=====================================================================================

   {
      test_ = "Row-major reset() function";

      using blaze::reset;
      using blaze::isDefault;

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 2UL );

      reset( sm(0,1) );

      checkRows    ( sm  , 3UL );
      checkColumns ( sm  , 2UL );
      checkNonZeros( sm  , 2UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( !isDefault( sm(0,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 0 )\n( -2 0 )\n(  0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  0 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major reset
   //=====================================================================================

   {
      test_ = "Row-major Submatrix::reset()";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 2UL );

      reset( sm );

      checkRows    ( sm  , 3UL );
      checkColumns ( sm  , 2UL );
      checkNonZeros( sm  , 0UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 7UL );

      if( !isDefault( sm ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  0 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "(  0  0 -3  0 )\n"
                                     "(  0  0  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major row-wise reset
   //=====================================================================================

   {
      test_ = "Row-major Submatrix::reset( size_t )";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 2UL );

      // Resetting the 0th row
      {
         reset( sm, 0UL );

         checkRows    ( sm  , 3UL );
         checkColumns ( sm  , 2UL );
         checkNonZeros( sm  , 2UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 9UL );

         if( sm(0,0) !=  0 || sm(0,1) != 0 ||
             sm(1,0) != -2 || sm(1,1) != 0 ||
             sm(2,0) !=  0 || sm(2,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 0th row failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0 0 )\n( -2 0 )\n(  0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st row
      {
         reset( sm, 1UL );

         checkRows    ( sm  , 3UL );
         checkColumns ( sm  , 2UL );
         checkNonZeros( sm  , 1UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 8UL );

         if( sm(0,0) != 0 || sm(0,1) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 0 ||
             sm(2,0) != 0 || sm(2,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st row failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd row
      {
         reset( sm, 2UL );

         checkRows    ( sm  , 3UL );
         checkColumns ( sm  , 2UL );
         checkNonZeros( sm  , 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 7UL );

         if( sm(0,0) != 0 || sm(0,1) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 0 ||
             sm(2,0) != 0 || sm(2,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd row failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major single element reset
   //=====================================================================================

   {
      test_ = "Column-major reset() function";

      using blaze::reset;
      using blaze::isDefault;

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 3UL );

      reset( sm(1,0) );

      checkRows    ( sm   , 2UL );
      checkColumns ( sm   , 3UL );
      checkNonZeros( sm   , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( !isDefault( sm(1,0) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -2 0 )\n( 0  0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 0 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  0  0  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major reset
   //=====================================================================================

   {
      test_ = "Column-major Submatrix::reset()";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 3UL );

      reset( sm );

      checkRows    ( sm   , 2UL );
      checkColumns ( sm   , 3UL );
      checkNonZeros( sm   , 0UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( !isDefault( sm ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 0 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  0  0  0 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major row-wise reset
   //=====================================================================================

   {
      test_ = "Column-major Submatrix::reset( size_t )";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 3UL );

      // Resetting the 0th column
      {
         reset( sm, 0UL );

         checkRows    ( sm   , 2UL );
         checkColumns ( sm   , 3UL );
         checkNonZeros( sm   , 2UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 9UL );

         if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 0 ||
             sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st column
      {
         reset( sm, 1UL );

         checkRows    ( sm   , 2UL );
         checkColumns ( sm   , 3UL );
         checkNonZeros( sm   , 1UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 8UL );

         if( sm(0,0) != 0 || sm(0,1) != 0 || sm(0,2) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st column failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd column
      {
         reset( sm, 2UL );

         checkRows    ( sm   , 2UL );
         checkColumns ( sm   , 3UL );
         checkNonZeros( sm   , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( sm(0,0) != 0 || sm(0,1) != 0 || sm(0,2) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd column failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() function with the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() function with the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testClear()
{
   //=====================================================================================
   // Row-major single element clear
   //=====================================================================================

   {
      test_ = "Row-major clear() function";

      using blaze::clear;
      using blaze::isDefault;

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 2UL );

      clear( sm(0,1) );

      checkRows    ( sm  , 3UL );
      checkColumns ( sm  , 2UL );
      checkNonZeros( sm  , 2UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( !isDefault( sm(0,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 0 )\n( -2 0 )\n(  0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  0 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major single element clear
   //=====================================================================================

   {
      test_ = "Column-major clear() function";

      using blaze::clear;
      using blaze::isDefault;

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 2UL, 3UL );

      clear( sm(1,0) );

      checkRows    ( sm   , 2UL );
      checkColumns ( sm   , 3UL );
      checkNonZeros( sm   , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( !isDefault( sm(1,0) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -2 0 )\n( 0  0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 0 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  0  0  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the Submatrix
// class template. Additionally, it performs a test of self-transpose via the \c trans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testTranspose()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via transpose()";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 3UL );

      transpose( sm );

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  5UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 0 ||
          sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 4 ||
          sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -2 0 )\n( 1  0 4 )\n( 0 -3 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -2 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 1 || mat_(2,1) !=  0 || mat_(2,2) != 4 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != -3 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 -2  0  0 )\n"
                                     "(  1  0  4  0 )\n"
                                     "(  0 -3  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via trans()";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 3UL );

      sm = trans( sm );

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  5UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 0 ||
          sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 4 ||
          sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -2 0 )\n( 1  0 4 )\n( 0 -3 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -2 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 1 || mat_(2,1) !=  0 || mat_(2,2) != 4 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != -3 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 -2  0  0 )\n"
                                     "(  1  0  4  0 )\n"
                                     "(  0 -3  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 3UL );

      transpose( sm );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  5UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
          sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -3 ||
          sm(2,0) !=  0 || sm(2,1) != 4 || sm(2,2) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  0  4  5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 1 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -2 || tmat_(1,2) != 0 || tmat_(1,3) != -3 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 4 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  1  0  7 )\n"
                                     "(  0 -2  0 -3 -8 )\n"
                                     "(  0  0  4  5  9 )\n"
                                     "(  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via trans()";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 3UL );

      sm = trans( sm );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  5UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
          sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -3 ||
          sm(2,0) !=  0 || sm(2,1) != 4 || sm(2,2) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  0  4  5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 1 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -2 || tmat_(1,2) != 0 || tmat_(1,3) != -3 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 4 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  1  0  7 )\n"
                                     "(  0 -2  0 -3 -8 )\n"
                                     "(  0  0  4  5  9 )\n"
                                     "(  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() member function of the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the Submatrix
// specialization. Additionally, it performs a test of self-transpose via the \c ctrans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testCTranspose()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via ctranspose()";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 3UL );

      ctranspose( sm );

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  5UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 0 ||
          sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 4 ||
          sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -2 0 )\n( 1  0 4 )\n( 0 -3 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -2 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 1 || mat_(2,1) !=  0 || mat_(2,2) != 4 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != -3 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 -2  0  0 )\n"
                                     "(  1  0  4  0 )\n"
                                     "(  0 -3  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via ctrans()";

      initialize();

      SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 3UL );

      sm = ctrans( sm );

      checkRows    ( sm  ,  3UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  5UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 0 ||
          sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 4 ||
          sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 -2 0 )\n( 1  0 4 )\n( 0 -3 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -2 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
          mat_(2,0) != 1 || mat_(2,1) !=  0 || mat_(2,2) != 4 || mat_(2,3) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != -3 || mat_(3,2) != 5 || mat_(3,3) != -6 ||
          mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 -2  0  0 )\n"
                                     "(  1  0  4  0 )\n"
                                     "(  0 -3  5 -6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 3UL );

      ctranspose( sm );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  5UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
          sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -3 ||
          sm(2,0) !=  0 || sm(2,1) != 4 || sm(2,2) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  0  4  5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 1 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -2 || tmat_(1,2) != 0 || tmat_(1,3) != -3 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 4 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  1  0  7 )\n"
                                     "(  0 -2  0 -3 -8 )\n"
                                     "(  0  0  4  5  9 )\n"
                                     "(  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via ctrans()";

      initialize();

      OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 3UL );

      sm = ctrans( sm );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  3UL );
      checkNonZeros( sm   ,  5UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
          sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -3 ||
          sm(2,0) !=  0 || sm(2,1) != 4 || sm(2,2) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  0  4  5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 1 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -2 || tmat_(1,2) != 0 || tmat_(1,3) != -3 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 4 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  1  0  7 )\n"
                                     "(  0 -2  0 -3 -8 )\n"
                                     "(  0  0  4  5  9 )\n"
                                     "(  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testIsDefault()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      using blaze::isDefault;

      initialize();

      // isDefault with default submatrix
      {
         SMT sm = blaze::submatrix( mat_, 0UL, 0UL, 1UL, 4UL );

         if( isDefault( sm(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix element: " << sm(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default submatrix
      {
         SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 1UL, 4UL );

         if( isDefault( sm(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix element: " << sm(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDefault() function";

      using blaze::isDefault;

      initialize();

      // isDefault with default submatrix
      {
         OSMT sm = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 1UL );

         if( isDefault( sm(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix element: " << sm(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default submatrix
      {
         OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 1UL );

         if( isDefault( sm(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix element: " << sm(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSame() function with the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSame() function with the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testIsSame()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major isSame() function";

      // isSame with matrix and matching submatrix
      {
         SMT sm = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 4UL );

         if( blaze::isSame( sm, mat_ ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( mat_, sm ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matrix and non-matching submatrix (different number of rows)
      {
         SMT sm = blaze::submatrix( mat_, 0UL, 0UL, 4UL, 4UL );

         if( blaze::isSame( sm, mat_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( mat_, sm ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matrix and non-matching submatrix (different number of columns)
      {
         SMT sm = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 3UL );

         if( blaze::isSame( sm, mat_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( mat_, sm ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matrix and non-matching submatrix (different row index)
      {
         SMT sm = blaze::submatrix( mat_, 1UL, 0UL, 4UL, 4UL );

         if( blaze::isSame( sm, mat_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( mat_, sm ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matrix and non-matching submatrix (different column index)
      {
         SMT sm = blaze::submatrix( mat_, 0UL, 1UL, 5UL, 3UL );

         if( blaze::isSame( sm, mat_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( mat_, sm ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching submatrices
      {
         SMT sm1 = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 4UL );
         SMT sm2 = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 4UL );

         if( blaze::isSame( sm1, sm2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching submatrices (different number of rows)
      {
         SMT sm1 = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 4UL );
         SMT sm2 = blaze::submatrix( mat_, 0UL, 0UL, 4UL, 4UL );

         if( blaze::isSame( sm1, sm2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching submatrices (different number of columns)
      {
         SMT sm1 = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 4UL );
         SMT sm2 = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 3UL );

         if( blaze::isSame( sm1, sm2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching submatrices (different row index)
      {
         SMT sm1 = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 4UL );
         SMT sm2 = blaze::submatrix( mat_, 1UL, 0UL, 4UL, 4UL );

         if( blaze::isSame( sm1, sm2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching submatrices (different column index)
      {
         SMT sm1 = blaze::submatrix( mat_, 0UL, 0UL, 5UL, 4UL );
         SMT sm2 = blaze::submatrix( mat_, 0UL, 1UL, 5UL, 3UL );

         if( blaze::isSame( sm1, sm2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major isSame() function";

      // isSame with matrix and matching submatrix
      {
         OSMT sm = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 5UL );

         if( blaze::isSame( sm, tmat_ ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( tmat_, sm ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matrix and non-matching submatrix (different number of rows)
      {
         OSMT sm = blaze::submatrix( tmat_, 0UL, 0UL, 3UL, 5UL );

         if( blaze::isSame( sm, tmat_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( tmat_, sm ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matrix and non-matching submatrix (different number of columns)
      {
         OSMT sm = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 4UL );

         if( blaze::isSame( sm, tmat_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( tmat_, sm ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matrix and non-matching submatrix (different row index)
      {
         OSMT sm = blaze::submatrix( tmat_, 1UL, 0UL, 3UL, 5UL );

         if( blaze::isSame( sm, tmat_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( tmat_, sm ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matrix and non-matching submatrix (different column index)
      {
         OSMT sm = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 4UL );

         if( blaze::isSame( sm, tmat_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( tmat_, sm ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << tmat_ << "\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching submatrices
      {
         OSMT sm1 = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 5UL );
         OSMT sm2 = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 5UL );

         if( blaze::isSame( sm1, sm2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching submatrices (different number of rows)
      {
         OSMT sm1 = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 5UL );
         OSMT sm2 = blaze::submatrix( tmat_, 0UL, 0UL, 3UL, 5UL );

         if( blaze::isSame( sm1, sm2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching submatrices (different number of columns)
      {
         OSMT sm1 = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 5UL );
         OSMT sm2 = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 4UL );

         if( blaze::isSame( sm1, sm2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching submatrices (different row index)
      {
         OSMT sm1 = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 5UL );
         OSMT sm2 = blaze::submatrix( tmat_, 1UL, 0UL, 3UL, 5UL );

         if( blaze::isSame( sm1, sm2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching submatrices (different column index)
      {
         OSMT sm1 = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 5UL );
         OSMT sm2 = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 4UL );

         if( blaze::isSame( sm1, sm2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First submatrix:\n" << sm1 << "\n"
                << "   Second submatrix:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c submatrix() function with the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c submatrix() function with the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      initialize();

      {
         SMT sm1 = blaze::submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = blaze::submatrix( sm1 , 1UL, 1UL, 3UL, 2UL );

         if( sm2(1,1) != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator access failed\n"
                << " Details:\n"
                << "   Result: " << sm2(1,1) << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }

         if( *sm2.begin(1UL) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *sm2.begin(1UL) << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         SMT sm1 = blaze::submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = blaze::submatrix( sm1 , 4UL, 1UL, 3UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm1 = blaze::submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = blaze::submatrix( sm1 , 1UL, 3UL, 3UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm1 = blaze::submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = blaze::submatrix( sm1 , 1UL, 1UL, 4UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm1 = blaze::submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = blaze::submatrix( sm1 , 1UL, 1UL, 3UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      initialize();

      {
         OSMT sm1 = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         OSMT sm2 = blaze::submatrix( sm1  , 1UL, 1UL, 2UL, 3UL );

         if( sm2(1,1) != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator access failed\n"
                << " Details:\n"
                << "   Result: " << sm2(1,1) << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }

         if( *sm2.begin(1UL) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *sm2.begin(1UL) << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         OSMT sm1 = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         OSMT sm2 = blaze::submatrix( sm1  , 3UL, 1UL, 2UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         OSMT sm1 = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         OSMT sm2 = blaze::submatrix( sm1  , 1UL, 4UL, 2UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         OSMT sm1 = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         OSMT sm2 = blaze::submatrix( sm1  , 1UL, 1UL, 3UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         OSMT sm1 = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         OSMT sm2 = blaze::submatrix( sm1  , 1UL, 1UL, 2UL, 4UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c row() function with the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c row() function with the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testRow()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      initialize();

      typedef blaze::Row<SMT>  RowType;

      SMT sm1 = blaze::submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
      RowType row1 = row( sm1, 1UL );

      if( row1[1] != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }

      if( *row1.begin() != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *row1.begin() << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      initialize();

      typedef blaze::Row<OSMT>  RowType;

      OSMT sm1 = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
      RowType row1 = row( sm1, 1UL );

      if( row1[1] != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }

      if( *row1.begin() != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *row1.begin() << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c column() function with the Submatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c column() function with the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testColumn()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      initialize();

      typedef blaze::Column<SMT>  ColumnType;

      SMT sm1 = blaze::submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
      ColumnType col1 = column( sm1, 1UL );

      if( col1[1] != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }

      if( *col1.begin() != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *col1.begin() << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      initialize();

      typedef blaze::Column<OSMT>  ColumnType;

      OSMT sm1 = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
      ColumnType col1 = column( sm1, 1UL );

      if( col1[1] != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }

      if( *col1.begin() != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *col1.begin() << "\n"
             << "   Expected result: 0\n";
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
void DenseUnalignedTest::initialize()
{
   // Initializing the row-major dynamic matrix
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

   // Initializing the column-major dynamic matrix
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

} // namespace submatrix

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
   std::cout << "   Running Submatrix dense unaligned test..." << std::endl;

   try
   {
      RUN_SUBMATRIX_DENSEUNALIGNED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Submatrix dense unaligned test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
