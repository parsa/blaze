//=================================================================================================
/*!
//  \file src/mathtest/sparsesubmatrix/ClassTest.cpp
//  \brief Source file for the SparseSubmatrix class test
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Views.h>
#include <blazetest/mathtest/sparsesubmatrix/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace sparsesubmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SparseSubmatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
   : mat_ ( 5UL, 4UL )
   , tmat_( 4UL, 5UL )
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testDivAssign();
   testFunctionCall();
   testIterator();
   testNonZeros();
   testReset();
   testAppend();
   testInsert();
   testErase();
   testReserve();
   testTrim();
   testScale();
   testFind();
   testLowerBound();
   testUpperBound();
   testIsDefault();
   testIsNan();
   testIsDiagonal();
   testIsSymmetric();
   testMinimum();
   testMaximum();
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
/*!\brief Test of the SparseSubmatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the SparseSubmatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix constructor";

      initialize();

      for( size_t row=0UL; row<mat_.rows(); ++row ) {
         for( size_t column=0UL; column<mat_.columns(); ++column ) {
            for( size_t m=1UL; (row+m)<mat_.rows(); ++m ) {
               for( size_t n=1UL; (column+n)<mat_.columns(); ++n )
               {
                  SMT sm = submatrix( mat_, row, column, m, n );

                  for( size_t i=0UL; i<m; ++i ) {
                     for( size_t j=0UL; j<n; ++j )
                     {
                        if( sm(i,j) != mat_(row+i,column+j) ) {
                           std::ostringstream oss;
                           oss << " Test: " << test_ << "\n"
                               << " Error: Setup of sparse submatrix failed\n"
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
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix constructor";

      initialize();

      for( size_t column=0UL; column<tmat_.columns(); ++column ) {
         for( size_t row=0UL; row<tmat_.rows(); ++row ) {
            for( size_t n=1UL; (column+n)<tmat_.columns(); ++n ) {
               for( size_t m=1UL; (row+m)<tmat_.rows(); ++m )
               {
                  TSMT sm = submatrix( tmat_, row, column, m, n );

                  for( size_t j=0UL; j<n; ++j ) {
                     for( size_t i=0UL; i<m; ++i )
                     {
                        if( sm(i,j) != tmat_(row+i,column+j) ) {
                           std::ostringstream oss;
                           oss << " Test: " << test_ << "\n"
                               << " Error: Setup of sparse submatrix failed\n"
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
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseSubmatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the SparseSubmatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix copy assignment (no aliasing)";

      initialize();

      MT mat( 5UL, 4UL, 3UL );
      mat(1,0) = 11;
      mat(2,0) = 12;
      mat(2,2) = 13;

      SMT sm = submatrix( mat, 1UL, 0UL, 2UL, 3UL );
      sm = submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

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
      test_ = "Row-major SparseSubmatrix copy assignment (aliasing)";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );
      sm = submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      /*
      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );
      */

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
      test_ = "Row-major/row-major dense matrix assignment";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
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
      test_ = "Row-major/column-major dense matrix assignment";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );
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

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 4UL );
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

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 4UL );
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
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix copy assignment (no aliasing)";

      initialize();

      TMT mat( 4UL, 5UL, 3UL );
      mat(0,1) = 11;
      mat(0,2) = 12;
      mat(2,2) = 13;

      TSMT sm = submatrix( mat, 0UL, 1UL, 3UL, 2UL );
      sm = submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

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
      test_ = "Column-major SparseSubmatrix copy assignment (aliasing)";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );
      sm = submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

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
      test_ = "Column-major/row-major dense matrix assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 0 );
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
      test_ = "Column-major/column-major dense matrix assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 4UL );
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 4UL );
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
/*!\brief Test of the SparseSubmatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the SparseSubmatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Row-major SparseSubmatrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix addition assignment (no aliasing)";

      initialize();

      MT mat( 5UL, 4UL, 3UL );
      mat(1,0) = 11;
      mat(2,0) = 12;
      mat(2,2) = 13;

      SMT sm = submatrix( mat, 1UL, 0UL, 2UL, 3UL );
      sm += submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

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
      test_ = "Row-major SparseSubmatrix addition assignment (aliasing)";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );
      sm += submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

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
      test_ = "Row-major/row-major dense matrix addition assignment";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
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
      test_ = "Row-major/column-major dense matrix addition assignment";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );
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

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 4UL );
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

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 4UL );
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
   // Column-major SparseSubmatrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix addition assignment (no aliasing)";

      initialize();

      TMT mat( 4UL, 5UL, 3UL );
      mat(0,1) = 11;
      mat(0,2) = 12;
      mat(2,2) = 13;

      TSMT sm = submatrix( mat, 0UL, 1UL, 3UL, 2UL );
      sm += submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

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
      test_ = "Column-major SparseSubmatrix addition assignment (aliasing)";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );
      sm += submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

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
      test_ = "Column-major/row-major dense matrix addition assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 0 );
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
      test_ = "Column-major/column-major dense matrix addition assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 4UL );
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 4UL );
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
/*!\brief Test of the SparseSubmatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the SparseSubmatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Row-major SparseSubmatrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix subtraction assignment (no aliasing)";

      initialize();

      MT mat( 5UL, 4UL, 3UL );
      mat(1,0) = 11;
      mat(2,0) = 12;
      mat(2,2) = 13;

      SMT sm = submatrix( mat, 1UL, 0UL, 2UL, 3UL );
      sm -= submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

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
      test_ = "Row-major SparseSubmatrix subtraction assignment (aliasing)";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );
      sm -= submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

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
      test_ = "Row-major/row-major dense matrix subtraction assignment";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
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
      test_ = "Row-major/column-major dense matrix subtraction assignment";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );
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

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 4UL );
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

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 4UL );
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
   // Column-major SparseSubmatrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix subtraction assignment (no aliasing)";

      initialize();

      TMT mat( 4UL, 5UL, 3UL );
      mat(0,1) = 11;
      mat(0,2) = 12;
      mat(2,2) = 13;

      TSMT sm = submatrix( mat, 0UL, 1UL, 3UL, 2UL );
      sm -= submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

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
      test_ = "Column-major SparseSubmatrix subtraction assignment (aliasing)";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );
      sm -= submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

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
      test_ = "Column-major/row-major dense matrix subtraction assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 0 );
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
      test_ = "Column-major/column-major dense matrix subtraction assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 4UL );
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 4UL );
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
/*!\brief Test of the SparseSubmatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the SparseSubmatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Row-major SparseSubmatrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix multiplication assignment (no aliasing)";

      initialize();

      MT mat( 5UL, 4UL, 4UL );
      mat(1,0) = 1;
      mat(1,1) = 1;
      mat(2,0) = 1;
      mat(2,1) = 1;

      SMT sm = submatrix( mat, 1UL, 0UL, 2UL, 2UL );
      sm *= submatrix( mat_, 2UL, 1UL, 2UL, 2UL );

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
      test_ = "Row-major SparseSubmatrix multiplication assignment (aliasing)";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 2UL );
      sm *= submatrix( mat_, 2UL, 1UL, 2UL, 2UL );

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
      test_ = "Row-major/row-major dense matrix multiplication assignment";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
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
      test_ = "Row-major/column-major dense matrix multiplication assignment";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
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

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
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

      SMT sm = submatrix( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
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
   // Row-major scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major scalar multiplication assignment";

      initialize();

      SMT sm = submatrix( mat_, 2UL, 0UL, 2UL, 3UL );

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
             << " Error: Multiplication assignment failed\n"
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
             << " Error: Multiplication assignment failed\n"
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
      test_ = "Row-major scalar multiplication assignment";

      initialize();

      SMT sm = submatrix( mat_, 2UL, 0UL, 3UL, 2UL );

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
             << " Error: Multiplication assignment failed\n"
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
             << " Error: Multiplication assignment failed\n"
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
   // Column-major SparseSubmatrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix multiplication assignment (no aliasing)";

      initialize();

      TMT mat( 4UL, 5UL, 4UL );
      mat(0,1) = 1;
      mat(0,2) = 1;
      mat(1,1) = 1;
      mat(1,2) = 1;

      TSMT sm = submatrix( mat, 0UL, 1UL, 2UL, 2UL );
      sm *= submatrix( tmat_, 1UL, 2UL, 2UL, 2UL );

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
      test_ = "Column-major SparseSubmatrix multiplication assignment (aliasing)";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );
      sm *= submatrix( tmat_, 1UL, 2UL, 2UL, 2UL );

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
      test_ = "Column-major/row-major dense matrix multiplication assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
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
      test_ = "Column-major/column-major dense matrix multiplication assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
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
   // Column-major scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major scalar multiplication assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 2UL, 3UL, 2UL );

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
             << " Error: Multiplication assignment failed\n"
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
             << " Error: Multiplication assignment failed\n"
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
      test_ = "Column-major scalar multiplication assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 2UL, 2UL, 3UL );

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
             << " Error: Multiplication assignment failed\n"
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
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  21 )\n"
                                     "( 0  1  0 12 -24 )\n"
                                     "( 0  0 -3  5   9 )\n"
                                     "( 0  0  0 -6  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseSubmatrix division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the SparseSubmatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testDivAssign()
{
   //=====================================================================================
   // Row-major scalar division assignment
   //=====================================================================================

   {
      test_ = "Row-major scalar division assignment";

      initialize();

      SMT sm = submatrix( mat_, 2UL, 0UL, 2UL, 3UL );

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
             << " Error: Division assignment failed\n"
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
             << " Error: Division assignment failed\n"
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
      test_ = "Row-major scalar division assignment";

      initialize();

      SMT sm = submatrix( mat_, 2UL, 0UL, 3UL, 2UL );

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
             << " Error: Division assignment failed\n"
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
             << " Error: Division assignment failed\n"
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
   // Column-major scalar division assignment
   //=====================================================================================

   {
      test_ = "Column-major scalar division assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 2UL, 3UL, 2UL );

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
             << " Error: Division assignment failed\n"
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
             << " Error: Division assignment failed\n"
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
      test_ = "Column-major scalar division assignment";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 2UL, 2UL, 3UL );

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
             << " Error: Division assignment failed\n"
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
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  14 )\n"
                                     "( 0  1  0  8 -16 )\n"
                                     "( 0  0 -3  5   9 )\n"
                                     "( 0  0  0 -6  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseSubmatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the SparseSubmatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::operator()";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 1UL, 3UL, 2UL );

      // Writing the first element
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

      // Writing the second element
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

      // Writing the third element
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
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::operator()";

      initialize();

      TSMT sm = submatrix( tmat_, 1UL, 1UL, 2UL, 3UL );

      // Writing the first element
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

      // Writing the second element
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

      // Writing the third element
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
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SparseSubmatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the SparseSubmatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 3UL, 3UL );

      // Counting the number of elements in 0th row
      {
         test_ = "Row-major iterator subtraction";

         const size_t number( sm.end(0) - sm.begin(0) );

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

      // Counting the number of elements in 1st row
      {
         test_ = "Row-major iterator subtraction";

         const size_t number( sm.end(1) - sm.begin(1) );

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

      // Counting the number of elements in 2nd row
      {
         test_ = "Row-major iterator subtraction";

         const size_t number( sm.end(2) - sm.begin(2) );

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

         SMT::ConstIterator it ( sm.cbegin(2) );
         SMT::ConstIterator end( sm.cend(2) );

         if( it == end || it->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it != sm.cend(2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Row-major assignment via Iterator";

         int value = 8;

         for( SMT::Iterator it=sm.begin(2); it!=sm.end(2); ++it ) {
            *it = value++;
         }

         if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
             sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -3 ||
             sm(2,0) !=  0 || sm(2,1) != 8 || sm(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  0  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  8 || mat_(3,2) !=  9 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  8  9 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         int value = 4;

         for( SMT::Iterator it=sm.begin(1); it!=sm.end(1); ++it ) {
            *it += value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 0 ||
             sm(1,0) != 2 || sm(1,1) != 0 || sm(1,2) != 2 ||
             sm(2,0) != 0 || sm(2,1) != 8 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 2 0 2 )\n( 0 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != 0 || mat_(0,3) !=  0 ||
             mat_(1,0) != 0 || mat_(1,1) !=  1 || mat_(1,2) != 0 || mat_(1,3) !=  0 ||
             mat_(2,0) != 2 || mat_(2,1) !=  0 || mat_(2,2) != 2 || mat_(2,3) !=  0 ||
             mat_(3,0) != 0 || mat_(3,1) !=  8 || mat_(3,2) != 9 || mat_(3,3) != -6 ||
             mat_(4,0) != 7 || mat_(4,1) != -8 || mat_(4,2) != 9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "(  2  0  2  0 )\n"
                                        "(  0  8  9 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         int value = 4;

         for( SMT::Iterator it=sm.begin(1); it!=sm.end(1); ++it ) {
            *it -= value++;
         }

         if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
             sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -3 ||
             sm(2,0) !=  0 || sm(2,1) != 8 || sm(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  0  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  8 || mat_(3,2) !=  9 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -3  0 )\n"
                                        "(  0  8  9 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         int value = 1;

         for( SMT::Iterator it=sm.begin(1); it!=sm.end(1); ++it ) {
            *it *= value++;
         }

         if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
             sm(1,0) != -2 || sm(1,1) != 0 || sm(1,2) != -6 ||
             sm(2,0) !=  0 || sm(2,1) != 8 || sm(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -6 )\n(  0  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -6 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  8 || mat_(3,2) !=  9 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -2  0 -6  0 )\n"
                                        "(  0  8  9 -6 )\n"
                                        "(  7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         for( SMT::Iterator it=sm.begin(1); it!=sm.end(1); ++it ) {
            *it /= 2;
         }

         if( sm(0,0) !=  0 || sm(0,1) != 1 || sm(0,2) !=  0 ||
             sm(1,0) != -1 || sm(1,1) != 0 || sm(1,2) != -3 ||
             sm(2,0) !=  0 || sm(2,1) != 8 || sm(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -1  0 -3 )\n(  0  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
             mat_(2,0) != -1 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) !=  8 || mat_(3,2) !=  9 || mat_(3,3) != -6 ||
             mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n(  0  0  0  0 )\n"
                                        "(  0  1  0  0 )\n"
                                        "( -1  0 -3  0 )\n"
                                        "(  0  8  9 -6 )\n"
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

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 3UL, 3UL );

      // Counting the number of elements in 0th column
      {
         test_ = "Column-major iterator subtraction";

         const size_t number( sm.end(0) - sm.begin(0) );

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

      // Counting the number of elements in 1st row
      {
         test_ = "Column-major iterator subtraction";

         const size_t number( sm.end(1) - sm.begin(1) );

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

      // Counting the number of elements in 2nd row
      {
         test_ = "Column-major iterator subtraction";

         const size_t number( sm.end(2) - sm.begin(2) );

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

         TSMT::ConstIterator it ( sm.cbegin(2) );
         TSMT::ConstIterator end( sm.cend(2) );

         if( it == end || it->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it != sm.cend(2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Column-major assignment via Iterator";

         int value = 8;

         for( TSMT::Iterator it=sm.begin(2); it!=sm.end(2); ++it ) {
            *it = value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 0 ||
             sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
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

         for( TSMT::Iterator it=sm.begin(1); it!=sm.end(1); ++it ) {
            *it += value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != 2 || sm(0,2) != 0 ||
             sm(1,0) != 1 || sm(1,1) != 0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != 2 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 2 0 )\n( 1 0 8 )\n( 0 2 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) != 0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 2 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0  2  0  7 )\n"
                                        "( 0  1  0  8 -8 )\n"
                                        "( 0  0  2  9  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         int value = 4;

         for( TSMT::Iterator it=sm.begin(1); it!=sm.end(1); ++it ) {
            *it -= value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 0 ||
             sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  1  0  8 -8 )\n"
                                        "( 0  0 -3  9  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         int value = 1;

         for( TSMT::Iterator it=sm.begin(1); it!=sm.end(1); ++it ) {
            *it *= value++;
         }

         if( sm(0,0) != 0 || sm(0,1) != -2 || sm(0,2) != 0 ||
             sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != -6 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 1  0  8 )\n( 0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -6 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  1  0  8 -8 )\n"
                                        "( 0  0 -6  9  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         for( TSMT::Iterator it=sm.begin(1); it!=sm.end(1); ++it ) {
            *it /= 2;
         }

         if( sm(0,0) != 0 || sm(0,1) != -1 || sm(0,2) != 0 ||
             sm(1,0) != 1 || sm(1,1) !=  0 || sm(1,2) != 8 ||
             sm(2,0) != 0 || sm(2,1) != -3 || sm(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 -1  0 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -1 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -1  0  7 )\n"
                                        "( 0  1  0  8 -8 )\n"
                                        "( 0  0 -3  9  9 )\n"
                                        "( 0  0  0 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the nonZeros member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of SparseSubmatrix. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::nonZeros()";

      initialize();

      // Initialization check
      SMT sm = submatrix( mat_, 1UL, 1UL, 2UL, 3UL );

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
      test_ = "Column-major SparseSubmatrix::nonZeros()";

      initialize();

      // Initialization check
      TSMT sm = submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );

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
/*!\brief Test of the reset member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   //=====================================================================================
   // Row-major reset
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::reset()";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 3UL, 2UL );

      sm.reset();

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
             << " Error: Subscript operator failed\n"
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
      test_ = "Row-major SparseSubmatrix::reset( size_t )";

      initialize();

      SMT sm = submatrix( mat_, 1UL, 0UL, 3UL, 2UL );

      // Resetting the 0th row
      {
         sm.reset( 0UL );

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
         sm.reset( 1UL );

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
         sm.reset( 2UL );

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
   // Column-major reset
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::reset()";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 2UL, 3UL );

      sm.reset();

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
             << " Error: Subscript operator failed\n"
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
      test_ = "Column-major SparseSubmatrix::reset( size_t )";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 2UL, 3UL );

      // Resetting the 0th column
      {
         sm.reset( 0UL );

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
         sm.reset( 1UL );

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
         sm.reset( 2UL );

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
/*!\brief Test of the append member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the append member function of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAppend()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::append()";

      // Appending with pre-allocation in each row
      {
         mat_.reset();

         // Initialization check
         SMT sm = submatrix( mat_, 0UL, 0UL, 4UL, 4UL );
         sm.reserve( 0UL, 2UL );
         sm.reserve( 2UL, 1UL );
         sm.reserve( 3UL, 2UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 0UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 0UL );
         checkNonZeros( sm, 3UL, 0UL );

         // Appending one non-zero element
         sm.append( 2UL, 1UL, 1 );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 1UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 1UL );
         checkNonZeros( sm, 3UL, 0UL );

         if( sm(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sm.append( 0UL, 0UL, 2 );
         sm.append( 0UL, 3UL, 3 );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 3UL );
         checkNonZeros( sm, 0UL, 2UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 1UL );
         checkNonZeros( sm, 3UL, 0UL );

         if( sm(2,1) != 1 || sm(0,0) != 2 || sm(0,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sm.append( 3UL, 1UL, 4 );
         sm.append( 3UL, 2UL, 5 );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 5UL );
         checkNonZeros( sm, 0UL, 2UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 1UL );
         checkNonZeros( sm, 3UL, 2UL );

         if( sm(2,1) != 1 || sm(0,0) != 2 || sm(0,3) != 3 ||
             sm(3,1) != 4 || sm(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 4 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         mat_.reset();

         // Initialization check
         SMT sm = submatrix( mat_, 0UL, 0UL, 4UL, 4UL );
         sm.reserve( 0UL, 2UL );
         sm.reserve( 2UL, 1UL );
         sm.reserve( 3UL, 2UL );

         // Appending one non-zero element
         sm.append( 0UL, 1UL, 1 );
         sm.finalize( 0UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 1UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 0UL );
         checkNonZeros( sm, 3UL, 0UL );

         if( sm(0,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sm.append( 1UL, 1UL, 2 );
         sm.append( 1UL, 3UL, 3 );
         sm.finalize( 1UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 3UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 2UL );
         checkNonZeros( sm, 2UL, 0UL );
         checkNonZeros( sm, 3UL, 0UL );

         if( sm(0,1) != 1 || sm(1,1) != 2 || sm(1,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 0 2 0 3 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sm.append( 3UL, 0UL, 4 );
         sm.append( 3UL, 1UL, 5 );
         sm.finalize( 1UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 5UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 2UL );
         checkNonZeros( sm, 2UL, 0UL );
         checkNonZeros( sm, 3UL, 2UL );

         if( sm(0,1) != 1 || sm(1,1) != 2 || sm(1,3) != 3 ||
             sm(3,0) != 4 || sm(3,1) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 0 2 0 3 )\n( 0 0 0 0 )\n( 4 5 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::append()";

      // Appending with pre-allocation in each row
      {
         tmat_.reset();

         // Initialization check
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 4UL, 4UL );
         sm.reserve( 0UL, 2UL );
         sm.reserve( 2UL, 1UL );
         sm.reserve( 3UL, 2UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 0UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 0UL );
         checkNonZeros( sm, 3UL, 0UL );

         // Appending one non-zero element
         sm.append( 1UL, 2UL, 1 );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 1UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 1UL );
         checkNonZeros( sm, 3UL, 0UL );

         if( sm(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sm.append( 0UL, 0UL, 2 );
         sm.append( 3UL, 0UL, 3 );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 3UL );
         checkNonZeros( sm, 0UL, 2UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 1UL );
         checkNonZeros( sm, 3UL, 0UL );

         if( sm(1,2) != 1 || sm(0,0) != 2 || sm(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sm.append( 1UL, 3UL, 4 );
         sm.append( 2UL, 3UL, 5 );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 5UL );
         checkNonZeros( sm, 0UL, 2UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 1UL );
         checkNonZeros( sm, 3UL, 2UL );

         if( sm(1,2) != 1 || sm(0,0) != 2 || sm(3,0) != 3 ||
             sm(1,3) != 4 || sm(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 1 4 )\n( 0 0 0 5 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         tmat_.reset();

         // Initialization check
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 4UL, 4UL );
         sm.reserve( 0UL, 2UL );
         sm.reserve( 2UL, 1UL );
         sm.reserve( 3UL, 2UL );

         // Appending one non-zero element
         sm.append( 1UL, 0UL, 1 );
         sm.finalize( 0UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 1UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 0UL );
         checkNonZeros( sm, 2UL, 0UL );
         checkNonZeros( sm, 3UL, 0UL );

         if( sm(1,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sm.append( 1UL, 1UL, 2 );
         sm.append( 3UL, 1UL, 3 );
         sm.finalize( 1UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 3UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 2UL );
         checkNonZeros( sm, 2UL, 0UL );
         checkNonZeros( sm, 3UL, 0UL );

         if( sm(1,0) != 1 || sm(1,1) != 2 || sm(3,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 2 0 0 )\n( 0 0 0 0 )\n( 0 3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sm.append( 0UL, 3UL, 4 );
         sm.append( 1UL, 3UL, 5 );
         sm.finalize( 1UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkCapacity( sm, 5UL );
         checkNonZeros( sm, 5UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 2UL );
         checkNonZeros( sm, 2UL, 0UL );
         checkNonZeros( sm, 3UL, 2UL );

         if( sm(1,0) != 1 || sm(1,1) != 2 || sm(3,1) != 3 ||
             sm(0,3) != 4 || sm(1,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 4 )\n( 1 2 0 5 )\n( 0 0 0 0 )\n( 0 3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the insert member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the insert member function of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testInsert()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::insert()";

      initialize();

      SMT sm = submatrix( mat_, 0UL, 1UL, 2UL, 3UL );

      // Inserting a non-zero element at the end of the 0th row
      sm.insert( 0UL, 2UL, 1 );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  2UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 11UL );

      if( sm(0,0) != 0 || sm(0,1) != 0 || sm(0,2) != 1 ||
          sm(1,0) != 1 || sm(1,1) != 0 || sm(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 1 )\n( 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Inserting a non-zero element at the beginning of the 0th row
      sm.insert( 0UL, 0UL, 2 );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  3UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( sm(0,0) != 2 || sm(0,1) != 0 || sm(0,2) != 1 ||
          sm(1,0) != 1 || sm(1,1) != 0 || sm(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 2 0 1 )\n( 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Inserting a non-zero element at the center of the 0th row
      sm.insert( 0UL, 1UL, 3 );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkNonZeros( sm  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 13UL );

      if( sm(0,0) != 2 || sm(0,1) != 3 || sm(0,2) != 1 ||
          sm(1,0) != 1 || sm(1,1) != 0 || sm(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 2 3 1 )\n( 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to insert an already existing element
      try {
         sm.insert( 1UL, 0UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 2 3 1 )\n( 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::insert()";

      initialize();

      TSMT sm = submatrix( tmat_, 1UL, 0UL, 3UL, 2UL );

      // Inserting a non-zero element at the end of the 0th column
      sm.insert( 2UL, 0UL, 1 );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 ||
          sm(1,0) != 0 || sm(1,1) != 0 ||
          sm(2,0) != 1 || sm(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Inserting a non-zero element at the beginning of the 0th column
      sm.insert( 0UL, 0UL, 2 );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  3UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( sm(0,0) != 2 || sm(0,1) != 1 ||
          sm(1,0) != 0 || sm(1,1) != 0 ||
          sm(2,0) != 1 || sm(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 2 1 )\n( 0 0 )\n( 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Inserting a non-zero element at the center of the 0th column
      sm.insert( 1UL, 0UL, 3 );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkNonZeros( sm   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 13UL );

      if( sm(0,0) != 2 || sm(0,1) != 1 ||
          sm(1,0) != 3 || sm(1,1) != 0 ||
          sm(2,0) != 1 || sm(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 2 1 )\n( 3 0 )\n( 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to insert an already existing element
      try {
         sm.insert( 0UL, 1UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 2 4 )\n( 3 0 )\n( 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the erase member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the erase member function of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::erase( size_t, size_t )";

      initialize();

      SMT sm = submatrix( mat_, 3UL, 1UL, 2UL, 3UL );

      // Erasing the non-zero element at the end of the 1st row
      sm.erase( 1UL, 2UL );

      checkRows    ( sm  , 2UL );
      checkColumns ( sm  , 3UL );
      checkNonZeros( sm  , 5UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( sm(0,0) !=  4 || sm(0,1) != 5 || sm(0,2) != -6 ||
          sm(1,0) != -8 || sm(1,1) != 9 || sm(1,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  4  5 -6 )\n( -8  9  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the 1st row
      sm.erase( 1UL, 0UL );

      checkRows    ( sm  , 2UL );
      checkColumns ( sm  , 3UL );
      checkNonZeros( sm  , 4UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 8UL );

      if( sm(0,0) != 4 || sm(0,1) != 5 || sm(0,2) != -6 ||
          sm(1,0) != 0 || sm(1,1) != 9 || sm(1,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 4  5 -6 )\n( 0  9  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the 1st row
      sm.erase( 1UL, 1UL );

      checkRows    ( sm  , 2UL );
      checkColumns ( sm  , 3UL );
      checkNonZeros( sm  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 7UL );

      if( sm(0,0) != 4 || sm(0,1) != 5 || sm(0,2) != -6 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 4  5 -6 )\n( 0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      sm.erase( 1UL, 2UL );

      checkRows    ( sm  , 2UL );
      checkColumns ( sm  , 3UL );
      checkNonZeros( sm  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 7UL );

      if( sm(0,0) != 4 || sm(0,1) != 5 || sm(0,2) != -6 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 4  5 -6 )\n( 0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::erase( size_t, Iterator )";

      initialize();

      SMT sm = submatrix( mat_, 3UL, 1UL, 2UL, 3UL );

      // Erasing the non-zero element at the end of the 1st row
      {
         SMT::Iterator pos = sm.erase( 1UL, sm.find( 1UL, 2UL ) );

         checkRows    ( sm  , 2UL );
         checkColumns ( sm  , 3UL );
         checkNonZeros( sm  , 5UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 9UL );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) !=  4 || sm(0,1) != 5 || sm(0,2) != -6 ||
             sm(1,0) != -8 || sm(1,1) != 9 || sm(1,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  4  5 -6 )\n( -8  9  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the 1st row
      {
         SMT::Iterator pos = sm.erase( 1UL, sm.find( 1UL, 0UL ) );

         checkRows    ( sm  , 2UL );
         checkColumns ( sm  , 3UL );
         checkNonZeros( sm  , 4UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 8UL );

         if( pos->value() != 9 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 9\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 4 || sm(0,1) != 5 || sm(0,2) != -6 ||
             sm(1,0) != 0 || sm(1,1) != 9 || sm(1,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 4  5 -6 )\n( 0  9  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the 1st row
      {
         SMT::Iterator pos = sm.erase( 1UL, sm.find( 1UL, 1UL ) );

         checkRows    ( sm  , 2UL );
         checkColumns ( sm  , 3UL );
         checkNonZeros( sm  , 3UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 7UL );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 4 || sm(0,1) != 5 || sm(0,2) != -6 ||
             sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 4  5 -6 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         SMT::Iterator pos = sm.erase( 1UL, sm.find( 1UL, 2UL ) );

         checkRows    ( sm  , 2UL );
         checkColumns ( sm  , 3UL );
         checkNonZeros( sm  , 3UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 7UL );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 4 || sm(0,1) != 5 || sm(0,2) != -6 ||
             sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 4  5 -6 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::erase( size_t, Iterator, Iterator )";

      initialize();

      SMT sm = submatrix( mat_, 3UL, 0UL, 2UL, 4UL );

      // Erasing the 0th row
      {
         SMT::Iterator pos = sm.erase( 0UL, sm.begin( 0UL ), sm.end( 0UL ) );

         checkRows    ( sm  , 2UL );
         checkColumns ( sm  , 4UL );
         checkNonZeros( sm  , 4UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 7UL );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 0 || sm(0,1) !=  0 || sm(0,2) != 0 || sm(0,3) !=  0 ||
             sm(1,0) != 7 || sm(1,1) != -8 || sm(1,2) != 9 || sm(1,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th row failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0  0  0  0 )\n( 7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 1st row
      {
         SMT::Iterator pos = sm.erase( 1UL, sm.begin( 1UL ), sm.find( 1UL, 2UL ) );

         checkRows    ( sm  , 2UL );
         checkColumns ( sm  , 4UL );
         checkNonZeros( sm  , 2UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 5UL );

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

         if( sm(0,0) != 0 || sm(0,1) != 0 || sm(0,2) != 0 || sm(0,3) !=  0 ||
             sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 9 || sm(1,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the first half of the 1st row failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0  0  0  0 )\n( 0  0  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 1st row
      {
         SMT::Iterator pos = sm.erase( 1UL, sm.find( 1UL, 2UL ), sm.end( 1UL ) );

         checkRows    ( sm  , 2UL );
         checkColumns ( sm  , 4UL );
         checkNonZeros( sm  , 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 3UL );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 0 || sm(0,1) != 0 || sm(0,2) != 0 || sm(0,3) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 0 || sm(1,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the second half of the 1st row failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         SMT::Iterator pos = sm.erase( 1UL, sm.begin( 1UL ), sm.begin( 1UL ) );

         checkRows    ( sm  , 2UL );
         checkColumns ( sm  , 4UL );
         checkNonZeros( sm  , 0UL );
         checkRows    ( mat_, 5UL );
         checkColumns ( mat_, 4UL );
         checkNonZeros( mat_, 3UL );

         if( pos != sm.begin( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 0 || sm(0,1) != 0 || sm(0,2) != 0 || sm(0,3) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 0 || sm(1,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::erase( size_t, size_t )";

      initialize();

      TSMT sm = submatrix( tmat_, 1UL, 3UL, 3UL, 2UL );

      // Erasing the non-zero element at the end of the 1st column
      sm.erase( 2UL, 1UL );

      checkRows    ( sm   , 3UL );
      checkColumns ( sm   , 2UL );
      checkNonZeros( sm   , 5UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( sm(0,0) !=  4 || sm(0,1) != -8 ||
          sm(1,0) !=  5 || sm(1,1) !=  9 ||
          sm(2,0) != -6 || sm(2,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  4 -8 )\n(  5  9 )\n( -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the 1st column
      sm.erase( 0UL, 1UL );

      checkRows    ( sm   , 3UL );
      checkColumns ( sm   , 2UL );
      checkNonZeros( sm   , 4UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 8UL );

      if( sm(0,0) !=  4 || sm(0,1) != 0 ||
          sm(1,0) !=  5 || sm(1,1) != 9 ||
          sm(2,0) != -6 || sm(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  4 0 )\n(  5 9 )\n( -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the 1st column
      sm.erase( 1UL, 1UL );

      checkRows    ( sm   , 3UL );
      checkColumns ( sm   , 2UL );
      checkNonZeros( sm   , 3UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( sm(0,0) !=  4 || sm(0,1) != 0 ||
          sm(1,0) !=  5 || sm(1,1) != 0 ||
          sm(2,0) != -6 || sm(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  4 0 )\n(  5 0 )\n( -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      sm.erase( 2UL, 1UL );

      checkRows    ( sm   , 3UL );
      checkColumns ( sm   , 2UL );
      checkNonZeros( sm   , 3UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( sm(0,0) !=  4 || sm(0,1) != 0 ||
          sm(1,0) !=  5 || sm(1,1) != 0 ||
          sm(2,0) != -6 || sm(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  4 0 )\n(  5 0 )\n( -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::erase( size_t, Iterator )";

      initialize();

      TSMT sm = submatrix( tmat_, 1UL, 3UL, 3UL, 2UL );

      // Erasing the non-zero element at the end of the 1st column
      {
         TSMT::Iterator pos = sm.erase( 1UL, sm.find( 2UL, 1UL ) );

         checkRows    ( sm   , 3UL );
         checkColumns ( sm   , 2UL );
         checkNonZeros( sm   , 5UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 9UL );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) !=  4 || sm(0,1) != -8 ||
             sm(1,0) !=  5 || sm(1,1) !=  9 ||
             sm(2,0) != -6 || sm(2,1) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  4 -8 )\n(  5  9 )\n( -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the 1st column
      {
         TSMT::Iterator pos = sm.erase( 1UL, sm.find( 0UL, 1UL ) );

         checkRows    ( sm   , 3UL );
         checkColumns ( sm   , 2UL );
         checkNonZeros( sm   , 4UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 8UL );

         if( pos->value() != 9 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 9\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) !=  4 || sm(0,1) != 0 ||
             sm(1,0) !=  5 || sm(1,1) != 9 ||
             sm(2,0) != -6 || sm(2,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  4 0 )\n(  5 9 )\n( -6 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the 1st column
      {
         TSMT::Iterator pos = sm.erase( 1UL, sm.find( 1UL, 1UL ) );

         checkRows    ( sm   , 3UL );
         checkColumns ( sm   , 2UL );
         checkNonZeros( sm   , 3UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) !=  4 || sm(0,1) != 0 ||
             sm(1,0) !=  5 || sm(1,1) != 0 ||
             sm(2,0) != -6 || sm(2,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  4 0 )\n(  5 0 )\n( -6 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         TSMT::Iterator pos = sm.erase( 1UL, sm.find( 2UL, 1UL ) );

         checkRows    ( sm   , 3UL );
         checkColumns ( sm   , 2UL );
         checkNonZeros( sm   , 3UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) !=  4 || sm(0,1) != 0 ||
             sm(1,0) !=  5 || sm(1,1) != 0 ||
             sm(2,0) != -6 || sm(2,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  4 0 )\n(  5 0 )\n( -6 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::erase( size_t, Iterator, Iterator )";

      initialize();

      TSMT sm = submatrix( tmat_, 0UL, 3UL, 4UL, 2UL );

      // Erasing the 0th column
      {
         TSMT::Iterator pos = sm.erase( 0UL, sm.begin( 0UL ), sm.end( 0UL ) );

         checkRows    ( sm   , 4UL );
         checkColumns ( sm   , 2UL );
         checkNonZeros( sm   , 4UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 0 || sm(0,1) !=  7 ||
             sm(1,0) != 0 || sm(1,1) != -8 ||
             sm(2,0) != 0 || sm(2,1) !=  9 ||
             sm(3,0) != 0 || sm(3,1) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0  7 )\n( 0 -8 )\n( 0  9 )\n( 0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 1st column
      {
         TSMT::Iterator pos = sm.erase( 1UL, sm.begin( 1UL ), sm.find( 2UL, 1UL ) );

         checkRows    ( sm   , 4UL );
         checkColumns ( sm   , 2UL );
         checkNonZeros( sm   , 2UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 5UL );

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

         if( sm(0,0) != 0 || sm(0,1) !=  0 ||
             sm(1,0) != 0 || sm(1,1) !=  0 ||
             sm(2,0) != 0 || sm(2,1) !=  9 ||
             sm(3,0) != 0 || sm(3,1) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0  0 )\n( 0  0 )\n( 0  9 )\n( 0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 1st column
      {
         TSMT::Iterator pos = sm.erase( 1UL, sm.find( 2UL, 1UL ), sm.end( 1UL ) );

         checkRows    ( sm   , 4UL );
         checkColumns ( sm   , 2UL );
         checkNonZeros( sm   , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 3UL );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 0 || sm(0,1) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 0 ||
             sm(2,0) != 0 || sm(2,1) != 0 ||
             sm(3,0) != 0 || sm(3,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         TSMT::Iterator pos = sm.erase( 1UL, sm.begin( 1UL ), sm.begin( 1UL ) );

         checkRows    ( sm   , 4UL );
         checkColumns ( sm   , 2UL );
         checkNonZeros( sm   , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 3UL );

         if( pos != sm.begin( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm(0,0) != 0 || sm(0,1) != 0 ||
             sm(1,0) != 0 || sm(1,1) != 0 ||
             sm(2,0) != 0 || sm(2,1) != 0 ||
             sm(3,0) != 0 || sm(3,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reserve member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reserve member function of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReserve()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::reserve()";

      MT mat( 3UL, 20UL );

      SMT sm = submatrix( mat, 1UL, 0UL, 1UL, 20UL );

      // Increasing the capacity of the matrix
      sm.reserve( 10UL );

      checkRows    ( sm,  1UL );
      checkColumns ( sm, 20UL );
      checkCapacity( sm, 10UL );
      checkNonZeros( sm,  0UL );

      // Further increasing the capacity of the matrix
      sm.reserve( 20UL );

      checkRows    ( sm,  1UL );
      checkColumns ( sm, 20UL );
      checkCapacity( sm, 20UL );
      checkNonZeros( sm,  0UL );
   }

   {
      test_ = "Row-major SparseSubmatrix::reserve( size_t )";

      MT mat( 3UL, 20UL );

      SMT sm = submatrix( mat, 1UL, 0UL, 1UL, 20UL );

      // Increasing the capacity of the row
      sm.reserve( 0UL, 10UL );

      checkRows    ( sm,  1UL );
      checkColumns ( sm, 20UL );
      checkCapacity( sm, 10UL );
      checkNonZeros( sm,  0UL );

      // Further increasing the capacity of the row
      sm.reserve( 0UL, 15UL );

      checkRows    ( sm,  1UL );
      checkColumns ( sm, 20UL );
      checkCapacity( sm, 15UL );
      checkNonZeros( sm,  0UL );
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::reserve()";

      TMT mat( 3UL, 20UL );

      TSMT sm = submatrix( mat, 1UL, 0UL, 1UL, 20UL );

      // Increasing the capacity of the matrix
      sm.reserve( 10UL );

      checkRows    ( sm,  1UL );
      checkColumns ( sm, 20UL );
      checkCapacity( sm, 10UL );
      checkNonZeros( sm,  0UL );

      // Further increasing the capacity of the matrix
      sm.reserve( 20UL );

      checkRows    ( sm,  1UL );
      checkColumns ( sm, 20UL );
      checkCapacity( sm, 20UL );
      checkNonZeros( sm,  0UL );
   }

   {
      test_ = "Columnt-major SparseSubmatrix::reserve( size_t )";

      TMT mat( 20UL, 3UL );

      TSMT sm = submatrix( mat, 0UL, 1UL, 20UL, 1UL );

      // Increasing the capacity of the column
      sm.reserve( 0UL, 10UL );

      checkRows    ( sm, 20UL );
      checkColumns ( sm,  1UL );
      checkCapacity( sm, 10UL );
      checkNonZeros( sm,  0UL );

      // Further increasing the capacity of the column
      sm.reserve( 0UL, 15UL );

      checkRows    ( sm, 20UL );
      checkColumns ( sm,  1UL );
      checkCapacity( sm, 15UL );
      checkNonZeros( sm,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the trim member functions of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the trim member functions of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::trim()";

      initialize();

      SMT sm = submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      // Increasing the row capacity of the matrix
      sm.reserve( 0UL, 10UL );
      sm.reserve( 1UL, 20UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkCapacity( sm  , 30UL );
      checkCapacity( sm  ,  0UL, 10UL );
      checkCapacity( sm  ,  1UL, 20UL );
      checkCapacity( mat_, 30UL );
      checkCapacity( mat_,  2UL, 10UL );
      checkCapacity( mat_,  3UL, 20UL );

      // Trimming the matrix
      sm.trim();

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkCapacity( sm  , 30UL );
      checkCapacity( sm  ,  0UL, sm.nonZeros( 0UL ) );
      checkCapacity( sm  ,  1UL, sm.nonZeros( 1UL ) );
      checkCapacity( mat_, 30UL );
      checkCapacity( mat_,  2UL, mat_.nonZeros( 2UL ) );
      checkCapacity( mat_,  3UL, mat_.nonZeros( 3UL ) );
   }

   {
      test_ = "Row-major CompressedMatrix::trim( size_t )";

      initialize();

      SMT sm = submatrix( mat_, 2UL, 1UL, 2UL, 3UL );

      // Increasing the row capacity of the matrix
      sm.reserve( 0UL, 10UL );
      sm.reserve( 1UL, 20UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkCapacity( sm  , 30UL );
      checkCapacity( sm  ,  0UL, 10UL );
      checkCapacity( sm  ,  1UL, 20UL );
      checkCapacity( mat_, 30UL );
      checkCapacity( mat_,  2UL, 10UL );
      checkCapacity( mat_,  3UL, 20UL );

      // Trimming the 0th row
      sm.trim( 0UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkCapacity( sm  , 30UL );
      checkCapacity( sm  ,  0UL, sm.nonZeros( 0UL ) );
      checkCapacity( sm  ,  1UL, 30UL - sm.nonZeros( 0UL ) );
      checkCapacity( mat_, 30UL );
      checkCapacity( mat_,  2UL, mat_.nonZeros( 2UL ) );
      checkCapacity( mat_,  3UL, 30UL - mat_.nonZeros( 2UL ) );

      // Trimming the 1st row
      sm.trim( 1UL );

      checkRows    ( sm  ,  2UL );
      checkColumns ( sm  ,  3UL );
      checkCapacity( sm  , 30UL );
      checkCapacity( sm  ,  0UL, sm.nonZeros( 0UL ) );
      checkCapacity( sm  ,  1UL, sm.nonZeros( 1UL ) );
      checkCapacity( mat_, 30UL );
      checkCapacity( mat_,  2UL, mat_.nonZeros( 2UL ) );
      checkCapacity( mat_,  3UL, mat_.nonZeros( 3UL ) );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::trim()";

      initialize();

      TSMT sm = submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

      // Increasing the row capacity of the matrix
      sm.reserve( 0UL, 10UL );
      sm.reserve( 1UL, 20UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkCapacity( sm   , 30UL );
      checkCapacity( sm   ,  0UL, 10UL );
      checkCapacity( sm   ,  1UL, 20UL );
      checkCapacity( tmat_, 30UL );
      checkCapacity( tmat_,  2UL, 10UL );
      checkCapacity( tmat_,  3UL, 20UL );

      // Trimming the matrix
      sm.trim();

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkCapacity( sm   , 30UL );
      checkCapacity( sm   ,  0UL, sm.nonZeros( 0UL ) );
      checkCapacity( sm   ,  1UL, sm.nonZeros( 1UL ) );
      checkCapacity( tmat_, 30UL );
      checkCapacity( tmat_,  2UL, tmat_.nonZeros( 2UL ) );
      checkCapacity( tmat_,  3UL, tmat_.nonZeros( 3UL ) );
   }

   {
      test_ = "Column-major CompressedMatrix::trim( size_t )";

      initialize();

      TSMT sm = submatrix( tmat_, 1UL, 2UL, 3UL, 2UL );

      // Increasing the row capacity of the matrix
      sm.reserve( 0UL, 10UL );
      sm.reserve( 1UL, 20UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkCapacity( sm   , 30UL );
      checkCapacity( sm   ,  0UL, 10UL );
      checkCapacity( sm   ,  1UL, 20UL );
      checkCapacity( tmat_, 30UL );
      checkCapacity( tmat_,  2UL, 10UL );
      checkCapacity( tmat_,  3UL, 20UL );

      // Trimming the 0th row
      sm.trim( 0UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkCapacity( sm   , 30UL );
      checkCapacity( sm   ,  0UL, sm.nonZeros( 0UL ) );
      checkCapacity( sm   ,  1UL, 30UL - sm.nonZeros( 0UL ) );
      checkCapacity( tmat_, 30UL );
      checkCapacity( tmat_,  2UL, tmat_.nonZeros( 2UL ) );
      checkCapacity( tmat_,  3UL, 30UL - tmat_.nonZeros( 2UL ) );

      // Trimming the 1st row
      sm.trim( 1UL );

      checkRows    ( sm   ,  3UL );
      checkColumns ( sm   ,  2UL );
      checkCapacity( sm   , 30UL );
      checkCapacity( sm   ,  0UL, sm.nonZeros( 0UL ) );
      checkCapacity( sm   ,  1UL, sm.nonZeros( 1UL ) );
      checkCapacity( tmat_, 30UL );
      checkCapacity( tmat_,  2UL, tmat_.nonZeros( 2UL ) );
      checkCapacity( tmat_,  3UL, tmat_.nonZeros( 3UL ) );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of SparseSubmatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScale()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::scale()";

      initialize();

      // Initialization check
      SMT sm = submatrix( mat_, 2UL, 1UL, 2UL, 2UL );

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
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::scale()";

      initialize();

      // Initialization check
      TSMT sm = submatrix( tmat_, 1UL, 2UL, 2UL, 2UL );

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
/*!\brief Test of the find member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the find member function of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::find()";

      typedef SMT::ConstIterator  ConstIterator;

      initialize();

      SMT sm = submatrix( mat_, 1UL, 1UL, 3UL, 2UL );

      checkRows    ( sm, 3UL );
      checkColumns ( sm, 2UL );
      checkNonZeros( sm, 4UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 1UL );
      checkNonZeros( sm, 2UL, 2UL );

      // Searching for the first element
      {
         ConstIterator pos( sm.find( 0UL, 0UL ) );

         if( pos == sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
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
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( sm.find( 1UL, 1UL ) );

         if( pos == sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 1 || pos->value() != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( sm.find( 1UL, 0UL ) );

         if( pos != sm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::find()";

      typedef TSMT::ConstIterator  ConstIterator;

      initialize();

      TSMT sm = submatrix( tmat_, 1UL, 1UL, 2UL, 3UL );

      checkRows    ( sm, 2UL );
      checkColumns ( sm, 3UL );
      checkNonZeros( sm, 4UL );
      checkNonZeros( sm, 0UL, 1UL );
      checkNonZeros( sm, 1UL, 1UL );
      checkNonZeros( sm, 2UL, 2UL );

      // Searching for the first element
      {
         ConstIterator pos( sm.find( 0UL, 0UL ) );

         if( pos == sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
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
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( sm.find( 1UL, 2UL ) );

         if( pos == sm.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 1 || pos->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 5\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( sm.find( 1UL, 0UL ) );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the lowerBound member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the lowerBound member function of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::lowerBound()";

      typedef SMT::ConstIterator  ConstIterator;

      SMT sm = submatrix( mat_, 1UL, 0UL, 1UL, 4UL );

      checkRows    ( sm, 1UL );
      checkColumns ( sm, 4UL );
      checkNonZeros( sm, 1UL );
      checkNonZeros( sm, 0UL, 1UL );

      // Determining the lower bound for position (0,0)
      {
         ConstIterator pos( sm.lowerBound( 0UL, 0UL ) );

         if( pos == sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
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
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (0,1)
      {
         ConstIterator pos( sm.lowerBound( 0UL, 1UL ) );

         if( pos == sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current submatrix:\n" << sm << "\n";
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
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (0,2)
      {
         ConstIterator pos( sm.lowerBound( 0UL, 2UL ) );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,2)\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::lowerBound()";

      typedef TSMT::ConstIterator  ConstIterator;

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 4UL, 1UL );

      checkRows    ( sm, 4UL );
      checkColumns ( sm, 1UL );
      checkNonZeros( sm, 1UL );
      checkNonZeros( sm, 0UL, 1UL );

      // Determining the lower bound for position (0,0)
      {
         ConstIterator pos( sm.lowerBound( 0UL, 0UL ) );

         if( pos == sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
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
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,0)
      {
         ConstIterator pos( sm.lowerBound( 1UL, 0UL ) );

         if( pos == sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
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
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,0)
      {
         ConstIterator pos( sm.lowerBound( 2UL, 0UL ) );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the upperBound member function of SparseSubmatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the upperBound member function of SparseSubmatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major SparseSubmatrix::upperBound()";

      typedef SMT::ConstIterator  ConstIterator;

      SMT sm = submatrix( mat_, 1UL, 0UL, 1UL, 4UL );

      checkRows    ( sm, 1UL );
      checkColumns ( sm, 4UL );
      checkNonZeros( sm, 1UL );
      checkNonZeros( sm, 0UL, 1UL );

      // Determining the upper bound for position (0,0)
      {
         ConstIterator pos( sm.upperBound( 0UL, 0UL ) );

         if( pos == sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
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
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (0,1)
      {
         ConstIterator pos( sm.upperBound( 0UL, 1UL ) );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (0,2)
      {
         ConstIterator pos( sm.upperBound( 0UL, 2UL ) );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,2)\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatrix::upperBound()";

      typedef TSMT::ConstIterator  ConstIterator;

      TSMT sm = submatrix( tmat_, 0UL, 1UL, 4UL, 1UL );

      checkRows    ( sm, 4UL );
      checkColumns ( sm, 1UL );
      checkNonZeros( sm, 1UL );
      checkNonZeros( sm, 0UL, 1UL );

      // Determining the upper bound for position (0,0)
      {
         ConstIterator pos( sm.upperBound( 0UL, 0UL ) );

         if( pos == sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
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
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,0)
      {
         ConstIterator pos( sm.upperBound( 1UL, 0UL ) );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,0)
      {
         ConstIterator pos( sm.upperBound( 2UL, 0UL ) );

         if( pos != sm.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,0)\n"
                << "   Current submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDefault function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDefault function with the SparseSubmatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      initialize();

      // isDefault with default submatrix
      {
         SMT sm = submatrix( mat_, 0UL, 0UL, 1UL, 4UL );

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
         SMT sm = submatrix( mat_, 1UL, 0UL, 1UL, 4UL );

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

      initialize();

      // isDefault with default submatrix
      {
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 4UL, 1UL );

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
         TSMT sm = submatrix( tmat_, 0UL, 1UL, 4UL, 1UL );

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
/*!\brief Test of the isnan function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isnan function with the SparseSubmatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsNan()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major isnan() function";

      typedef blaze::CompressedMatrix<float,blaze::rowMajor>  MatrixType;
      typedef blaze::SparseSubmatrix<MatrixType>              SubmatrixType;

      MatrixType mat( mat_ );

      // isnan with empty 2x2 matrix
      {
         SubmatrixType sm = submatrix( mat, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );

         if( blaze::isnan( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with filled 2x3 matrix
      {
         SubmatrixType sm = submatrix( mat, 2UL, 1UL, 2UL, 3UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 4UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 3UL );

         if( blaze::isnan( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
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
      test_ = "Column-major isnan() function";

      typedef blaze::CompressedMatrix<float,blaze::columnMajor>  MatrixType;
      typedef blaze::SparseSubmatrix<MatrixType>                 SubmatrixType;

      MatrixType mat( tmat_ );

      // isnan with empty 2x2 matrix
      {
         SubmatrixType sm = submatrix( mat, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );

         if( blaze::isnan( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with filled 3x2 matrix
      {
         SubmatrixType sm = submatrix( mat, 1UL, 2UL, 3UL, 2UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 4UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 3UL );

         if( blaze::isnan( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDiagonal function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDiagonal function with the SparseSubmatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDiagonal()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDiagonal() function";

      initialize();
      mat_(0,0) = 11;
      mat_(2,0) =  0;

      // Non-quadratic submatrix
      {
         SMT sm = submatrix( mat_, 0UL, 0UL, 2UL, 3UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 2UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );

         if( isDiagonal( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         SMT sm = submatrix( mat_, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );

         if( isDiagonal( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         SMT sm = submatrix( mat_, 0UL, 0UL, 3UL, 3UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 3UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );
         checkNonZeros( sm, 2UL, 1UL );

         if( isDiagonal( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-diagonal matrix
      {
         SMT sm = submatrix( mat_, 0UL, 0UL, 4UL, 4UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkNonZeros( sm, 6UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );
         checkNonZeros( sm, 2UL, 1UL );
         checkNonZeros( sm, 3UL, 3UL );

         if( isDiagonal( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
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
      test_ = "Column-major isDiagonal() function";

      initialize();
      tmat_(0,0) = 11;
      tmat_(0,2) =  0;

      // Non-quadratic submatrix
      {
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 3UL, 2UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 2UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );

         if( isDiagonal( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         TSMT sm = submatrix( tmat_, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );

         if( isDiagonal( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 3UL, 3UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 3UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );
         checkNonZeros( sm, 2UL, 1UL );

         if( isDiagonal( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-diagonal matrix
      {
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 4UL, 4UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkNonZeros( sm, 6UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );
         checkNonZeros( sm, 2UL, 1UL );
         checkNonZeros( sm, 3UL, 3UL );

         if( isDiagonal( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isSymmetric function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isSymmetric function with the SparseSubmatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsSymmetric()
{
   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major isSymmetric() function";

      initialize();
      mat_(0,0) = 11;
      mat_(2,0) =  0;
      mat_(2,3) =  5;
      mat_(3,1) =  0;

      // Non-quadratic matrix
      {
         SMT sm = submatrix( mat_, 0UL, 0UL, 2UL, 3UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 2UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );

         if( isSymmetric( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         SMT sm = submatrix( mat_, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );

         if( isSymmetric( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         SMT sm = submatrix( mat_, 0UL, 0UL, 3UL, 3UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 3UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );
         checkNonZeros( sm, 2UL, 1UL );

         if( isSymmetric( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix
      {
         SMT sm = submatrix( mat_, 1UL, 0UL, 4UL, 4UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkNonZeros( sm, 9UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 2UL );
         checkNonZeros( sm, 2UL, 2UL );
         checkNonZeros( sm, 3UL, 4UL );

         if( isSymmetric( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Symmetric matrix
      {
         SMT sm = submatrix( mat_, 0UL, 0UL, 4UL, 4UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkNonZeros( sm, 6UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );
         checkNonZeros( sm, 2UL, 2UL );
         checkNonZeros( sm, 3UL, 2UL );

         if( isSymmetric( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
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
      test_ = "Column-major isSymmetric() function";

      initialize();
      tmat_(0,0) = 11;
      tmat_(0,2) =  0;
      tmat_(3,2) =  5;
      tmat_(1,3) =  0;

      // Non-quadratic matrix
      {
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 3UL, 2UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 2UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );

         if( isSymmetric( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         TSMT sm = submatrix( tmat_, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );
         checkNonZeros( sm, 0UL, 0UL );
         checkNonZeros( sm, 1UL, 0UL );

         if( isSymmetric( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 3UL, 3UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 3UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );
         checkNonZeros( sm, 2UL, 1UL );

         if( isSymmetric( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix
      {
         TSMT sm = submatrix( tmat_, 0UL, 1UL, 4UL, 4UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkNonZeros( sm, 9UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 2UL );
         checkNonZeros( sm, 2UL, 2UL );
         checkNonZeros( sm, 3UL, 4UL );

         if( isSymmetric( sm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Symmetric matrix
      {
         TSMT sm = submatrix( tmat_, 0UL, 0UL, 4UL, 4UL );

         checkRows    ( sm, 4UL );
         checkColumns ( sm, 4UL );
         checkNonZeros( sm, 6UL );
         checkNonZeros( sm, 0UL, 1UL );
         checkNonZeros( sm, 1UL, 1UL );
         checkNonZeros( sm, 2UL, 2UL );
         checkNonZeros( sm, 3UL, 2UL );

         if( isSymmetric( sm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the min function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the min function with the SparseSubmatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMinimum()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major min() function";

      initialize();

      // Attempt to find the minimum in an empty submatrix
      {
         SMT sm = submatrix( mat_, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );

         const int minimum = min( sm );

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

      // Attempt to find the minimum in a partially filled submatrix
      {
         SMT sm = submatrix( mat_, 1UL, 1UL, 2UL, 3UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 2UL );

         const int minimum = min( sm );

         if( minimum != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the minimum in a fully filled submatrix
      {
         SMT sm = submatrix( mat_, 3UL, 1UL, 2UL, 3UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 6UL );

         const int minimum = min( sm );

         if( minimum != -8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Third computation failed\n"
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

      // Attempt to find the minimum in an empty submatrix
      {
         TSMT sm = submatrix( tmat_, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );

         const int minimum = min( sm );

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

      // Attempt to find the minimum in a partially filled submatrix
      {
         TSMT sm = submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 2UL );

         const int minimum = min( sm );

         if( minimum != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the minimum in a fully filled submatrix
      {
         TSMT sm = submatrix( tmat_, 1UL, 3UL, 3UL, 2UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 6UL );

         const int minimum = min( sm );

         if( minimum != -8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Third computation failed\n"
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
/*!\brief Test of the max function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the max function with the SparseSubmatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMaximum()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major max() function";

      initialize();

      // Attempt to find the maximum in an empty submatrix
      {
         SMT sm = submatrix( mat_, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );

         const int maximum = max( sm );

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

      // Attempt to find the maximum in a partially filled submatrix
      {
         SMT sm = submatrix( mat_, 1UL, 1UL, 2UL, 3UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 2UL );

         const int maximum = max( sm );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the maximum in a fully filled submatrix
      {
         SMT sm = submatrix( mat_, 3UL, 1UL, 2UL, 3UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 3UL );
         checkNonZeros( sm, 6UL );

         const int maximum = max( sm );

         if( maximum != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Third computation failed\n"
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

      // Attempt to find the maximum in an empty submatrix
      {
         TSMT sm = submatrix( tmat_, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( sm, 2UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 0UL );

         const int maximum = max( sm );

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

      // Attempt to find the maximum in a partially filled submatrix
      {
         TSMT sm = submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 2UL );

         const int maximum = max( sm );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the maximum in a fully filled submatrix
      {
         TSMT sm = submatrix( tmat_, 1UL, 3UL, 3UL, 2UL );

         checkRows    ( sm, 3UL );
         checkColumns ( sm, 2UL );
         checkNonZeros( sm, 6UL );

         const int maximum = max( sm );

         if( maximum != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Third computation failed\n"
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
/*!\brief Test of the submatrix function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the submatrix function with the SparseSubmatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      initialize();

      {
         SMT sm1 = submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = submatrix( sm1 , 1UL, 1UL, 3UL, 2UL );

         if( sm2(1,1) != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator access failed\n"
                << " Details:\n"
                << "   Result: " << sm2(1,1) << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm2.begin(1UL)->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << sm2.begin(1UL)->value() << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         SMT sm1 = submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = submatrix( sm1 , 4UL, 1UL, 3UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm1 = submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = submatrix( sm1 , 1UL, 3UL, 3UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm1 = submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = submatrix( sm1 , 1UL, 1UL, 4UL, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         SMT sm1 = submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
         SMT sm2 = submatrix( sm1 , 1UL, 1UL, 3UL, 3UL );

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
         TSMT sm1 = submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         TSMT sm2 = submatrix( sm1  , 1UL, 1UL, 2UL, 3UL );

         if( sm2(1,1) != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator access failed\n"
                << " Details:\n"
                << "   Result: " << sm2(1,1) << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }

         if( sm2.begin(1UL)->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << sm2.begin(1UL)->value() << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         TSMT sm1 = submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         TSMT sm2 = submatrix( sm1  , 3UL, 1UL, 2UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         TSMT sm1 = submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         TSMT sm2 = submatrix( sm1  , 1UL, 4UL, 2UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         TSMT sm1 = submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         TSMT sm2 = submatrix( sm1  , 1UL, 1UL, 3UL, 3UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         TSMT sm1 = submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
         TSMT sm2 = submatrix( sm1  , 1UL, 1UL, 2UL, 4UL );

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
/*!\brief Test of the row function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the row function with the SparseSubmatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testRow()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      initialize();

      typedef blaze::SparseRow<SMT>  RowType;

      SMT sm1 = submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
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

      if( row1.begin()->value() != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << row1.begin()->value() << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      initialize();

      typedef blaze::SparseRow<TSMT>  RowType;

      TSMT sm1 = submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
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

      if( row1.begin()->value() != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << row1.begin()->value() << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the column function with the SparseSubmatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the column function with the SparseSubmatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testColumn()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      initialize();

      typedef blaze::SparseColumn<SMT>  ColumnType;

      SMT sm1 = submatrix( mat_, 1UL, 1UL, 4UL, 3UL );
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

      if( col1.begin()->value() != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << col1.begin()->value() << "\n"
             << "   Expected result: -3\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      initialize();

      typedef blaze::SparseColumn<TSMT>  ColumnType;

      TSMT sm1 = submatrix( tmat_, 1UL, 1UL, 3UL, 4UL );
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

      if( col1.begin()->value() != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << col1.begin()->value() << "\n"
             << "   Expected result: -3\n";
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

} // namespace sparsesubmatrix

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
   std::cout << "   Running SparseSubmatrix class test..." << std::endl;

   try
   {
      RUN_SPARSESUBMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SparseSubmatrix class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
