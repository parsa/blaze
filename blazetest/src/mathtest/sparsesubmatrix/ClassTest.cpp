//=================================================================================================
/*!
//  \file src/mathtest/sparsesubmatrix/ClassTest.cpp
//  \brief Source file for the SparseSubmatrix class test
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
*/
//=================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <blaze/math/DynamicMatrix.h>
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
                  SMT submatrix = sub( mat_, row, column, m, n );

                  for( size_t i=0UL; i<m; ++i ) {
                     for( size_t j=0UL; j<n; ++j )
                     {
                        if( submatrix(i,j) != mat_(row+i,column+j) ) {
                           std::ostringstream oss;
                           oss << " Test: " << test_ << "\n"
                               << " Error: Setup of sparse submatrix failed\n"
                               << " Details:\n"
                               << "   Index of first row    = " << row << "\n"
                               << "   Index of first column = " << column << "\n"
                               << "   Number of rows        = " << m << "\n"
                               << "   Number of columns     = " << n << "\n"
                               << "   Submatrix:\n" << submatrix << "\n"
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
                  TSMT submatrix = sub( tmat_, row, column, m, n );

                  for( size_t j=0UL; j<n; ++j ) {
                     for( size_t i=0UL; i<m; ++i )
                     {
                        if( submatrix(i,j) != tmat_(row+i,column+j) ) {
                           std::ostringstream oss;
                           oss << " Test: " << test_ << "\n"
                               << " Error: Setup of sparse submatrix failed\n"
                               << " Details:\n"
                               << "   Index of first row    = " << row << "\n"
                               << "   Index of first column = " << column << "\n"
                               << "   Number of rows        = " << m << "\n"
                               << "   Number of columns     = " << n << "\n"
                               << "   Submatrix:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat, 1UL, 0UL, 2UL, 3UL );
      submatrix = sub( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 10UL );
      checkRows    ( mat      ,  5UL );
      checkColumns ( mat      ,  4UL );
      checkNonZeros( mat      ,  4UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != -3 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 4 || submatrix(1,1) !=  5 || submatrix(1,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );
      submatrix = sub( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != -3 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 4 || submatrix(1,1) !=  5 || submatrix(1,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      submatrix = mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 11 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 || submatrix(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      submatrix = mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 11 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 || submatrix(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      submatrix = mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 11 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 || submatrix(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      submatrix = mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 11 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 || submatrix(1,2) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( mat, 0UL, 1UL, 3UL, 2UL );
      submatrix = sub( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 10UL );
      checkRows    ( mat      ,  4UL );
      checkColumns ( mat      ,  5UL );
      checkNonZeros( mat      ,  4UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) !=  4 ||
          submatrix(1,0) != -3 || submatrix(1,1) !=  5 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );
      submatrix = sub( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) !=  4 ||
          submatrix(1,0) != -3 || submatrix(1,1) !=  5 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 0 );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      submatrix = mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 ||
          submatrix(1,0) != 11 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      submatrix = mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 ||
          submatrix(1,0) != 11 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      submatrix = mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 ||
          submatrix(1,0) != 11 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      submatrix = mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 ||
          submatrix(1,0) != 11 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat, 1UL, 0UL, 2UL, 3UL );
      submatrix += sub( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  5UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 10UL );
      checkRows    ( mat      ,  5UL );
      checkColumns ( mat      ,  4UL );
      checkNonZeros( mat      ,  5UL );

      if( submatrix(0,0) != 11 || submatrix(0,1) != -3 || submatrix(0,2) != 0 ||
          submatrix(1,0) != 16 || submatrix(1,1) !=  5 || submatrix(1,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );
      submatrix += sub( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != -2 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 2 || submatrix(1,1) !=  5 || submatrix(1,2) != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      submatrix += mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 10 || submatrix(1,1) != 13 || submatrix(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      submatrix += mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 10 || submatrix(1,1) != 13 || submatrix(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      submatrix += mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 10 || submatrix(1,1) != 13 || submatrix(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = 11;
      mat(1,0) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      submatrix += mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 10 || submatrix(1,1) != 13 || submatrix(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( mat, 0UL, 1UL, 3UL, 2UL );
      submatrix += sub( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  5UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 10UL );
      checkRows    ( mat      ,  4UL );
      checkColumns ( mat      ,  5UL );
      checkNonZeros( mat      ,  5UL );

      if( submatrix(0,0) != 11 || submatrix(0,1) != 16 ||
          submatrix(1,0) != -3 || submatrix(1,1) !=  5 ||
          submatrix(2,0) !=  0 || submatrix(2,1) !=  7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );
      submatrix += sub( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) !=  2 ||
          submatrix(1,0) != -2 || submatrix(1,1) !=  5 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 0 );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      submatrix += mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 10 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      submatrix += mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 10 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      submatrix += mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 10 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 13;
      mat(2,1) = 14;

      submatrix += mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 10 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat, 1UL, 0UL, 2UL, 3UL );
      submatrix -= sub( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  5UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 10UL );
      checkRows    ( mat      ,  5UL );
      checkColumns ( mat      ,  4UL );
      checkNonZeros( mat      ,  5UL );

      if( submatrix(0,0) != 11 || submatrix(0,1) !=  3 || submatrix(0,2) !=  0 ||
          submatrix(1,0) !=  8 || submatrix(1,1) != -5 || submatrix(1,2) != 19 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );
      submatrix -= sub( mat_, 2UL, 1UL, 2UL, 3UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) !=  4 || submatrix(0,2) != 0 ||
          submatrix(1,0) != -6 || submatrix(1,1) != -5 || submatrix(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      submatrix -= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 10 || submatrix(1,1) != 13 || submatrix(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      submatrix -= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 10 || submatrix(1,1) != 13 || submatrix(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      submatrix -= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 10 || submatrix(1,1) != 13 || submatrix(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 3UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 4UL );
      mat(0,1) = -11;
      mat(1,0) = -12;
      mat(1,1) = -13;
      mat(1,2) = -14;

      submatrix -= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 12 || submatrix(0,2) !=  0 ||
          submatrix(1,0) != 10 || submatrix(1,1) != 13 || submatrix(1,2) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( mat, 0UL, 1UL, 3UL, 2UL );
      submatrix -= sub( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  5UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 10UL );
      checkRows    ( mat      ,  4UL );
      checkColumns ( mat      ,  5UL );
      checkNonZeros( mat      ,  5UL );

      if( submatrix(0,0) != 11 || submatrix(0,1) !=  8 ||
          submatrix(1,0) !=  3 || submatrix(1,1) != -5 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 19 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );
      submatrix -= sub( tmat_, 1UL, 2UL, 3UL, 2UL );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != -6 ||
          submatrix(1,0) != 4 || submatrix(1,1) != -5 ||
          submatrix(2,0) != 0 || submatrix(2,1) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 0 );
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      submatrix -= mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 10 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      submatrix -= mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 10 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      submatrix -= mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 10 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 4UL );
      mat(1,0) = -11;
      mat(0,1) = -12;
      mat(1,1) = -13;
      mat(2,1) = -14;

      submatrix -= mat;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 10 ||
          submatrix(1,0) != 12 || submatrix(1,1) != 13 ||
          submatrix(2,0) !=  0 || submatrix(2,1) != 11 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat, 1UL, 0UL, 2UL, 2UL );
      submatrix *= sub( mat_, 2UL, 1UL, 2UL, 2UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 10UL );
      checkRows    ( mat      ,  5UL );
      checkColumns ( mat      ,  4UL );
      checkNonZeros( mat      ,  4UL );

      if( submatrix(0,0) != 4 || submatrix(0,1) != 2 ||
          submatrix(1,0) != 4 || submatrix(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 2UL );
      submatrix *= sub( mat_, 2UL, 1UL, 2UL, 2UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  3UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) != 4 || submatrix(0,1) != 5 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      submatrix *= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 12UL );

      if( submatrix(0,0) != 13 || submatrix(0,1) != 14 ||
          submatrix(1,0) != 22 || submatrix(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      submatrix *= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 12UL );

      if( submatrix(0,0) != 13 || submatrix(0,1) != 14 ||
          submatrix(1,0) != 22 || submatrix(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      submatrix *= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 12UL );

      if( submatrix(0,0) != 13 || submatrix(0,1) != 14 ||
          submatrix(1,0) != 22 || submatrix(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = -11;
      mat(0,1) = -12;
      mat(1,0) =  13;
      mat(1,1) =  14;

      submatrix *= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 12UL );

      if( submatrix(0,0) != 13 || submatrix(0,1) != 14 ||
          submatrix(1,0) != 22 || submatrix(1,1) != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 2UL, 0UL, 2UL, 3UL );

      submatrix *= 3;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 10UL );

      if( submatrix(0,0) != -6 || submatrix(0,1) !=  0 || submatrix(0,2) != -9 ||
          submatrix(1,0) !=  0 || submatrix(1,1) != 12 || submatrix(1,2) != 15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 2UL, 0UL, 3UL, 2UL );

      submatrix *= 3;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 10UL );

      if( submatrix(0,0) != -6 || submatrix(0,1) !=   0 ||
          submatrix(1,0) !=  0 || submatrix(1,1) !=  12 ||
          submatrix(2,0) != 21 || submatrix(2,1) != -24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( mat, 0UL, 1UL, 2UL, 2UL );
      submatrix *= sub( tmat_, 1UL, 2UL, 2UL, 2UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 10UL );
      checkRows    ( mat      ,  4UL );
      checkColumns ( mat      ,  5UL );
      checkNonZeros( mat      ,  4UL );

      if( submatrix(0,0) != -3 || submatrix(0,1) != 9 ||
          submatrix(1,0) != -3 || submatrix(1,1) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 2UL, 2UL );
      submatrix *= sub( tmat_, 1UL, 2UL, 2UL, 2UL );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  3UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) != 6 || submatrix(0,1) != -10 ||
          submatrix(1,0) != 0 || submatrix(1,1) !=   4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      submatrix *= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 12UL );

      if( submatrix(0,0) != 26 || submatrix(0,1) != 28 ||
          submatrix(1,0) != 11 || submatrix(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      submatrix *= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 12UL );

      if( submatrix(0,0) != 26 || submatrix(0,1) != 28 ||
          submatrix(1,0) != 11 || submatrix(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      submatrix *= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 12UL );

      if( submatrix(0,0) != 26 || submatrix(0,1) != 28 ||
          submatrix(1,0) != 11 || submatrix(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 2UL, 2UL );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  11;
      mat(0,1) =  12;
      mat(1,0) = -13;
      mat(1,1) = -14;

      submatrix *= mat;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 12UL );

      if( submatrix(0,0) != 26 || submatrix(0,1) != 28 ||
          submatrix(1,0) != 11 || submatrix(1,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 2UL, 3UL, 2UL );

      submatrix *= 3;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 10UL );

      if( submatrix(0,0) != -6 || submatrix(0,1) !=  0 ||
          submatrix(1,0) !=  0 || submatrix(1,1) != 12 ||
          submatrix(2,0) != -9 || submatrix(2,1) != 15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 2UL, 2UL, 3UL );

      submatrix *= 3;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 10UL );

      if( submatrix(0,0) != -6 || submatrix(0,1) !=  0 || submatrix(0,2) !=  21 ||
          submatrix(1,0) !=  0 || submatrix(1,1) != 12 || submatrix(1,2) != -24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 2UL, 0UL, 2UL, 3UL );

      submatrix /= 0.5;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 10UL );

      if( submatrix(0,0) != -4 || submatrix(0,1) != 0 || submatrix(0,2) != -6 ||
          submatrix(1,0) !=  0 || submatrix(1,1) != 8 || submatrix(1,2) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 2UL, 0UL, 3UL, 2UL );

      submatrix /= 0.5;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 10UL );

      if( submatrix(0,0) != -4 || submatrix(0,1) !=   0 ||
          submatrix(1,0) !=  0 || submatrix(1,1) !=   8 ||
          submatrix(2,0) != 14 || submatrix(2,1) != -16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 2UL, 3UL, 2UL );

      submatrix /= 0.5;

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 10UL );

      if( submatrix(0,0) != -4 || submatrix(0,1) !=  0 ||
          submatrix(1,0) !=  0 || submatrix(1,1) !=  8 ||
          submatrix(2,0) != -6 || submatrix(2,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 2UL, 2UL, 3UL );

      submatrix /= 0.5;

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 10UL );

      if( submatrix(0,0) != -4 || submatrix(0,1) != 0 || submatrix(0,2) !=  14 ||
          submatrix(1,0) !=  0 || submatrix(1,1) != 8 || submatrix(1,2) != -16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 1UL, 3UL, 2UL );

      // Writing the first element
      {
         submatrix(1,0) = 9;

         checkRows    ( submatrix,  3UL );
         checkColumns ( submatrix,  2UL );
         checkNonZeros( submatrix,  5UL );
         checkNonZeros( submatrix,  0UL, 1UL );
         checkNonZeros( submatrix,  1UL, 2UL );
         checkNonZeros( submatrix,  2UL, 2UL );
         checkRows    ( mat_     ,  5UL );
         checkColumns ( mat_     ,  4UL );
         checkNonZeros( mat_     , 11UL );

         if( submatrix(0,0) != 1 || submatrix(0,1) !=  0 ||
             submatrix(1,0) != 9 || submatrix(1,1) != -3 ||
             submatrix(2,0) != 4 || submatrix(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
         submatrix(2,0) = 0;

         checkRows    ( submatrix,  3UL );
         checkColumns ( submatrix,  2UL );
         checkNonZeros( submatrix,  4UL );
         checkNonZeros( submatrix,  0UL, 1UL );
         checkNonZeros( submatrix,  1UL, 2UL );
         checkNonZeros( submatrix,  2UL, 1UL );
         checkRows    ( mat_     ,  5UL );
         checkColumns ( mat_     ,  4UL );
         checkNonZeros( mat_     , 10UL );

         if( submatrix(0,0) != 1 || submatrix(0,1) !=  0 ||
             submatrix(1,0) != 9 || submatrix(1,1) != -3 ||
             submatrix(2,0) != 0 || submatrix(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
         submatrix(1,1) = 11;

         checkRows    ( submatrix,  3UL );
         checkColumns ( submatrix,  2UL );
         checkNonZeros( submatrix,  4UL );
         checkNonZeros( submatrix,  0UL, 1UL );
         checkNonZeros( submatrix,  1UL, 2UL );
         checkNonZeros( submatrix,  2UL, 1UL );
         checkRows    ( mat_     ,  5UL );
         checkColumns ( mat_     ,  4UL );
         checkNonZeros( mat_     , 10UL );

         if( submatrix(0,0) != 1 || submatrix(0,1) !=  0 ||
             submatrix(1,0) != 9 || submatrix(1,1) != 11 ||
             submatrix(2,0) != 0 || submatrix(2,1) !=  5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 1UL, 1UL, 2UL, 3UL );

      // Writing the first element
      {
         submatrix(0,1) = 9;

         checkRows    ( submatrix,  2UL );
         checkColumns ( submatrix,  3UL );
         checkNonZeros( submatrix,  5UL );
         checkNonZeros( submatrix,  0UL, 1UL );
         checkNonZeros( submatrix,  1UL, 2UL );
         checkNonZeros( submatrix,  2UL, 2UL );
         checkRows    ( tmat_    ,  4UL );
         checkColumns ( tmat_    ,  5UL );
         checkNonZeros( tmat_    , 11UL );

         if( submatrix(0,0) != 1 || submatrix(0,1) !=  9 || submatrix(0,2) != 4 ||
             submatrix(1,0) != 0 || submatrix(1,1) != -3 || submatrix(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
         submatrix(0,2) = 0;

         checkRows    ( submatrix,  2UL );
         checkColumns ( submatrix,  3UL );
         checkNonZeros( submatrix,  4UL );
         checkNonZeros( submatrix,  0UL, 1UL );
         checkNonZeros( submatrix,  1UL, 2UL );
         checkNonZeros( submatrix,  2UL, 1UL );
         checkRows    ( tmat_    ,  4UL );
         checkColumns ( tmat_    ,  5UL );
         checkNonZeros( tmat_    , 10UL );

         if( submatrix(0,0) != 1 || submatrix(0,1) !=  9 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != -3 || submatrix(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
         submatrix(1,1) = 11;

         checkRows    ( submatrix,  2UL );
         checkColumns ( submatrix,  3UL );
         checkNonZeros( submatrix,  4UL );
         checkNonZeros( submatrix,  0UL, 1UL );
         checkNonZeros( submatrix,  1UL, 2UL );
         checkNonZeros( submatrix,  2UL, 1UL );
         checkRows    ( tmat_    ,  4UL );
         checkColumns ( tmat_    ,  5UL );
         checkNonZeros( tmat_    , 10UL );

         if( submatrix(0,0) != 1 || submatrix(0,1) !=  9 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 11 || submatrix(1,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 3UL, 3UL );

      // Counting the number of elements in 0th row
      {
         test_ = "Row-major iterator subtraction";

         const size_t number( submatrix.end(0) - submatrix.begin(0) );

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

         const size_t number( submatrix.end(1) - submatrix.begin(1) );

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

         const size_t number( submatrix.end(2) - submatrix.begin(2) );

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

         SMT::ConstIterator it( submatrix.cbegin(2) );

         if( it->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Unexpected iterator behavior\n"
                << " Details:\n"
                << "   Current value : " << it->value() << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Unexpected iterator behavior\n"
                << " Details:\n"
                << "   Current value : " << it->value() << "\n"
                << "   Expected value: 5\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it != submatrix.cend(2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator end\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Row-major assignment via Iterator";

         int value = 8;

         for( SMT::Iterator it=submatrix.begin(2); it!=submatrix.end(2); ++it ) {
            *it = value++;
         }

         if( submatrix(0,0) !=  0 || submatrix(0,1) != 1 || submatrix(0,2) !=  0 ||
             submatrix(1,0) != -2 || submatrix(1,1) != 0 || submatrix(1,2) != -3 ||
             submatrix(2,0) !=  0 || submatrix(2,1) != 8 || submatrix(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

         for( SMT::Iterator it=submatrix.begin(1); it!=submatrix.end(1); ++it ) {
            *it += value++;
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != 1 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 2 || submatrix(1,1) != 0 || submatrix(1,2) != 2 ||
             submatrix(2,0) != 0 || submatrix(2,1) != 8 || submatrix(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
                << " Error: Assignment via iterator failed\n"
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

         for( SMT::Iterator it=submatrix.begin(1); it!=submatrix.end(1); ++it ) {
            *it -= value++;
         }

         if( submatrix(0,0) !=  0 || submatrix(0,1) != 1 || submatrix(0,2) !=  0 ||
             submatrix(1,0) != -2 || submatrix(1,1) != 0 || submatrix(1,2) != -3 ||
             submatrix(2,0) !=  0 || submatrix(2,1) != 8 || submatrix(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         int value = 1;

         for( SMT::Iterator it=submatrix.begin(1); it!=submatrix.end(1); ++it ) {
            *it *= value++;
         }

         if( submatrix(0,0) !=  0 || submatrix(0,1) != 1 || submatrix(0,2) !=  0 ||
             submatrix(1,0) != -2 || submatrix(1,1) != 0 || submatrix(1,2) != -6 ||
             submatrix(2,0) !=  0 || submatrix(2,1) != 8 || submatrix(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
                << " Error: Assignment via iterator failed\n"
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

         for( SMT::Iterator it=submatrix.begin(1); it!=submatrix.end(1); ++it ) {
            *it /= 2;
         }

         if( submatrix(0,0) !=  0 || submatrix(0,1) != 1 || submatrix(0,2) !=  0 ||
             submatrix(1,0) != -1 || submatrix(1,1) != 0 || submatrix(1,2) != -3 ||
             submatrix(2,0) !=  0 || submatrix(2,1) != 8 || submatrix(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
                << " Error: Assignment via iterator failed\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 3UL, 3UL );

      // Counting the number of elements in 0th column
      {
         test_ = "Column-major iterator subtraction";

         const size_t number( submatrix.end(0) - submatrix.begin(0) );

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

         const size_t number( submatrix.end(1) - submatrix.begin(1) );

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

         const size_t number( submatrix.end(2) - submatrix.begin(2) );

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

         TSMT::ConstIterator it( submatrix.cbegin(2) );

         if( it->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Unexpected iterator behavior\n"
                << " Details:\n"
                << "   Current value : " << it->value() << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Unexpected iterator behavior\n"
                << " Details:\n"
                << "   Current value : " << it->value() << "\n"
                << "   Expected value: 5\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it != submatrix.cend(2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator end\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Column-major assignment via Iterator";

         int value = 8;

         for( TSMT::Iterator it=submatrix.begin(2); it!=submatrix.end(2); ++it ) {
            *it = value++;
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != -2 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 1 || submatrix(1,1) !=  0 || submatrix(1,2) != 8 ||
             submatrix(2,0) != 0 || submatrix(2,1) != -3 || submatrix(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

         for( TSMT::Iterator it=submatrix.begin(1); it!=submatrix.end(1); ++it ) {
            *it += value++;
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != 2 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 1 || submatrix(1,1) != 0 || submatrix(1,2) != 8 ||
             submatrix(2,0) != 0 || submatrix(2,1) != 2 || submatrix(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 2 0 )\n( 1 0 8 )\n( 0 2 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) != 0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != 2 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) != 0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
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

         for( TSMT::Iterator it=submatrix.begin(1); it!=submatrix.end(1); ++it ) {
            *it -= value++;
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != -2 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 1 || submatrix(1,1) !=  0 || submatrix(1,2) != 8 ||
             submatrix(2,0) != 0 || submatrix(2,1) != -3 || submatrix(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         int value = 1;

         for( TSMT::Iterator it=submatrix.begin(1); it!=submatrix.end(1); ++it ) {
            *it *= value++;
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != -2 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 1 || submatrix(1,1) !=  0 || submatrix(1,2) != 8 ||
             submatrix(2,0) != 0 || submatrix(2,1) != -6 || submatrix(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 1  0  8 )\n( 0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -6 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
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

         for( TSMT::Iterator it=submatrix.begin(1); it!=submatrix.end(1); ++it ) {
            *it /= 2;
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != -1 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 1 || submatrix(1,1) !=  0 || submatrix(1,2) != 8 ||
             submatrix(2,0) != 0 || submatrix(2,1) != -3 || submatrix(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 -1  0 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -1 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  8 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
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
      SMT submatrix = sub( mat_, 1UL, 1UL, 2UL, 3UL );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 2UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 1UL );

      if( submatrix(0,0) != 1 || submatrix(0,1) !=  0 || submatrix(0,2) != 0 ||
          submatrix(1,0) != 0 || submatrix(1,1) != -3 || submatrix(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 1  0 0 )\n( 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse submatrix
      submatrix(1,1) = 0;

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 1UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 0UL );

      if( submatrix(0,0) != 1 || submatrix(0,1) != 0 || submatrix(0,2) != 0 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      mat_(2,3) = 5;

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 2UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 1UL );

      if( submatrix(0,0) != 1 || submatrix(0,1) != 0 || submatrix(0,2) != 0 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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
      TSMT submatrix = sub( tmat_, 1UL, 1UL, 3UL, 2UL );

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 2UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 1UL );

      if( submatrix(0,0) != 1 || submatrix(0,1) !=  0 ||
          submatrix(1,0) != 0 || submatrix(1,1) != -3 ||
          submatrix(2,0) != 0 || submatrix(2,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 1  0 )\n( 0 -3 )\n( 0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse submatrix
      submatrix(1,1) = 0;

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 1UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 0UL );

      if( submatrix(0,0) != 1 || submatrix(0,1) != 0 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 0 ||
          submatrix(2,0) != 0 || submatrix(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      tmat_(3,2) = 5;

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 2UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 1UL );

      if( submatrix(0,0) != 1 || submatrix(0,1) != 0 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 0 ||
          submatrix(2,0) != 0 || submatrix(2,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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
      test_ = "Row-major SparseSubmatix::reset()";

      initialize();

      SMT submatrix = sub( mat_, 1UL, 0UL, 3UL, 2UL );

      submatrix.reset();

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 0UL );
      checkRows    ( mat_     , 5UL );
      checkColumns ( mat_     , 4UL );
      checkNonZeros( mat_     , 7UL );

      if( !isDefault( submatrix ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 3UL, 2UL );

      // Resetting the 0th row
      {
         submatrix.reset( 0UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 2UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 9UL );

         if( submatrix(0,0) !=  0 || submatrix(0,1) != 0 ||
             submatrix(1,0) != -2 || submatrix(1,1) != 0 ||
             submatrix(2,0) !=  0 || submatrix(2,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 0th row failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n(  0 0 )\n( -2 0 )\n(  0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st row
      {
         submatrix.reset( 1UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 1UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 8UL );

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 ||
             submatrix(2,0) != 0 || submatrix(2,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st row failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd row
      {
         submatrix.reset( 2UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 7UL );

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 ||
             submatrix(2,0) != 0 || submatrix(2,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd row failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major reset
   //=====================================================================================

   {
      test_ = "Column-major SparseSubmatix::reset()";

      initialize();

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 2UL, 3UL );

      submatrix.reset();

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 0UL );
      checkRows    ( tmat_    , 4UL );
      checkColumns ( tmat_    , 5UL );
      checkNonZeros( tmat_    , 7UL );

      if( !isDefault( submatrix ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 2UL, 3UL );

      // Resetting the 0th column
      {
         submatrix.reset( 0UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 2UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 9UL );

         if( submatrix(0,0) != 0 || submatrix(0,1) != -2 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) !=  0 || submatrix(1,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st column
      {
         submatrix.reset( 1UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 1UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 8UL );

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st column failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd column
      {
         submatrix.reset( 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 7UL );

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 || submatrix(0,2) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd column failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
         SMT submatrix = sub( mat_, 0UL, 0UL, 4UL, 4UL );
         submatrix.reserve( 0UL, 2UL );
         submatrix.reserve( 2UL, 1UL );
         submatrix.reserve( 3UL, 2UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 0UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 0UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         // Appending one non-zero element
         submatrix.append( 2UL, 1UL, 1 );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 1UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 1UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         if( submatrix(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         submatrix.append( 0UL, 0UL, 2 );
         submatrix.append( 0UL, 3UL, 3 );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL, 2UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 1UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         if( submatrix(2,1) != 1 || submatrix(0,0) != 2 || submatrix(0,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         submatrix.append( 3UL, 1UL, 4 );
         submatrix.append( 3UL, 2UL, 5 );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 5UL );
         checkNonZeros( submatrix, 0UL, 2UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 1UL );
         checkNonZeros( submatrix, 3UL, 2UL );

         if( submatrix(2,1) != 1 || submatrix(0,0) != 2 || submatrix(0,3) != 3 ||
             submatrix(3,1) != 4 || submatrix(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 4 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         mat_.reset();

         // Initialization check
         SMT submatrix = sub( mat_, 0UL, 0UL, 4UL, 4UL );
         submatrix.reserve( 0UL, 2UL );
         submatrix.reserve( 2UL, 1UL );
         submatrix.reserve( 3UL, 2UL );

         // Appending one non-zero element
         submatrix.append( 0UL, 1UL, 1 );
         submatrix.finalize( 0UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 1UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 0UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         if( submatrix(0,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         submatrix.append( 1UL, 1UL, 2 );
         submatrix.append( 1UL, 3UL, 3 );
         submatrix.finalize( 1UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 2UL );
         checkNonZeros( submatrix, 2UL, 0UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         if( submatrix(0,1) != 1 || submatrix(1,1) != 2 || submatrix(1,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 0 2 0 3 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         submatrix.append( 3UL, 0UL, 4 );
         submatrix.append( 3UL, 1UL, 5 );
         submatrix.finalize( 1UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 5UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 2UL );
         checkNonZeros( submatrix, 2UL, 0UL );
         checkNonZeros( submatrix, 3UL, 2UL );

         if( submatrix(0,1) != 1 || submatrix(1,1) != 2 || submatrix(1,3) != 3 ||
             submatrix(3,0) != 4 || submatrix(3,1) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 4UL, 4UL );
         submatrix.reserve( 0UL, 2UL );
         submatrix.reserve( 2UL, 1UL );
         submatrix.reserve( 3UL, 2UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 0UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 0UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         // Appending one non-zero element
         submatrix.append( 1UL, 2UL, 1 );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 1UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 1UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         if( submatrix(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         submatrix.append( 0UL, 0UL, 2 );
         submatrix.append( 3UL, 0UL, 3 );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL, 2UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 1UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         if( submatrix(1,2) != 1 || submatrix(0,0) != 2 || submatrix(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         submatrix.append( 1UL, 3UL, 4 );
         submatrix.append( 2UL, 3UL, 5 );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 5UL );
         checkNonZeros( submatrix, 0UL, 2UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 1UL );
         checkNonZeros( submatrix, 3UL, 2UL );

         if( submatrix(1,2) != 1 || submatrix(0,0) != 2 || submatrix(3,0) != 3 ||
             submatrix(1,3) != 4 || submatrix(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 1 4 )\n( 0 0 0 5 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         tmat_.reset();

         // Initialization check
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 4UL, 4UL );
         submatrix.reserve( 0UL, 2UL );
         submatrix.reserve( 2UL, 1UL );
         submatrix.reserve( 3UL, 2UL );

         // Appending one non-zero element
         submatrix.append( 1UL, 0UL, 1 );
         submatrix.finalize( 0UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 1UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 0UL );
         checkNonZeros( submatrix, 2UL, 0UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         if( submatrix(1,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         submatrix.append( 1UL, 1UL, 2 );
         submatrix.append( 3UL, 1UL, 3 );
         submatrix.finalize( 1UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 2UL );
         checkNonZeros( submatrix, 2UL, 0UL );
         checkNonZeros( submatrix, 3UL, 0UL );

         if( submatrix(1,0) != 1 || submatrix(1,1) != 2 || submatrix(3,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 2 0 0 )\n( 0 0 0 0 )\n( 0 3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         submatrix.append( 0UL, 3UL, 4 );
         submatrix.append( 1UL, 3UL, 5 );
         submatrix.finalize( 1UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkCapacity( submatrix, 5UL );
         checkNonZeros( submatrix, 5UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 2UL );
         checkNonZeros( submatrix, 2UL, 0UL );
         checkNonZeros( submatrix, 3UL, 2UL );

         if( submatrix(1,0) != 1 || submatrix(1,1) != 2 || submatrix(3,1) != 3 ||
             submatrix(0,3) != 4 || submatrix(1,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 0UL, 1UL, 2UL, 3UL );

      // Inserting a non-zero element at the end of the 0th row
      submatrix.insert( 0UL, 2UL, 1 );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  2UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 11UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != 0 || submatrix(0,2) != 1 ||
          submatrix(1,0) != 1 || submatrix(1,1) != 0 || submatrix(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 0 0 1 )\n( 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Inserting a non-zero element at the beginning of the 0th row
      submatrix.insert( 0UL, 0UL, 2 );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  3UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 12UL );

      if( submatrix(0,0) != 2 || submatrix(0,1) != 0 || submatrix(0,2) != 1 ||
          submatrix(1,0) != 1 || submatrix(1,1) != 0 || submatrix(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 2 0 1 )\n( 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Inserting a non-zero element at the center of the 0th row
      submatrix.insert( 0UL, 1UL, 3 );

      checkRows    ( submatrix,  2UL );
      checkColumns ( submatrix,  3UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( mat_     ,  5UL );
      checkColumns ( mat_     ,  4UL );
      checkNonZeros( mat_     , 13UL );

      if( submatrix(0,0) != 2 || submatrix(0,1) != 3 || submatrix(0,2) != 1 ||
          submatrix(1,0) != 1 || submatrix(1,1) != 0 || submatrix(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 2 3 1 )\n( 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to insert an already existing element
      try {
         submatrix.insert( 1UL, 0UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 1UL, 0UL, 3UL, 2UL );

      // Inserting a non-zero element at the end of the 0th column
      submatrix.insert( 2UL, 0UL, 1 );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  2UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 11UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != 1 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 0 ||
          submatrix(2,0) != 1 || submatrix(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Inserting a non-zero element at the beginning of the 0th column
      submatrix.insert( 0UL, 0UL, 2 );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  3UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 12UL );

      if( submatrix(0,0) != 2 || submatrix(0,1) != 1 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 0 ||
          submatrix(2,0) != 1 || submatrix(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 2 1 )\n( 0 0 )\n( 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Inserting a non-zero element at the center of the 0th column
      submatrix.insert( 1UL, 0UL, 3 );

      checkRows    ( submatrix,  3UL );
      checkColumns ( submatrix,  2UL );
      checkNonZeros( submatrix,  4UL );
      checkRows    ( tmat_    ,  4UL );
      checkColumns ( tmat_    ,  5UL );
      checkNonZeros( tmat_    , 13UL );

      if( submatrix(0,0) != 2 || submatrix(0,1) != 1 ||
          submatrix(1,0) != 3 || submatrix(1,1) != 0 ||
          submatrix(2,0) != 1 || submatrix(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 2 1 )\n( 3 0 )\n( 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to insert an already existing element
      try {
         submatrix.insert( 0UL, 1UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 3UL, 1UL, 2UL, 3UL );

      // Erasing the non-zero element at the end of the 1st row
      submatrix.erase( 1UL, 2UL );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 5UL );
      checkRows    ( mat_     , 5UL );
      checkColumns ( mat_     , 4UL );
      checkNonZeros( mat_     , 9UL );

      if( submatrix(0,0) !=  4 || submatrix(0,1) != 5 || submatrix(0,2) != -6 ||
          submatrix(1,0) != -8 || submatrix(1,1) != 9 || submatrix(1,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n(  4  5 -6 )\n( -8  9  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the 1st row
      submatrix.erase( 1UL, 0UL );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 4UL );
      checkRows    ( mat_     , 5UL );
      checkColumns ( mat_     , 4UL );
      checkNonZeros( mat_     , 8UL );

      if( submatrix(0,0) != 4 || submatrix(0,1) != 5 || submatrix(0,2) != -6 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 9 || submatrix(1,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 4  5 -6 )\n( 0  9  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the 1st row
      submatrix.erase( 1UL, 1UL );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 3UL );
      checkRows    ( mat_     , 5UL );
      checkColumns ( mat_     , 4UL );
      checkNonZeros( mat_     , 7UL );

      if( submatrix(0,0) != 4 || submatrix(0,1) != 5 || submatrix(0,2) != -6 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 4  5 -6 )\n( 0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      submatrix.erase( 1UL, 2UL );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 3UL );
      checkRows    ( mat_     , 5UL );
      checkColumns ( mat_     , 4UL );
      checkNonZeros( mat_     , 7UL );

      if( submatrix(0,0) != 4 || submatrix(0,1) != 5 || submatrix(0,2) != -6 ||
          submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 3UL, 1UL, 2UL, 3UL );

      // Erasing the non-zero element at the end of the 1st row
      {
         SMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 1UL, 2UL ) );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 5UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 9UL );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) !=  4 || submatrix(0,1) != 5 || submatrix(0,2) != -6 ||
             submatrix(1,0) != -8 || submatrix(1,1) != 9 || submatrix(1,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n(  4  5 -6 )\n( -8  9  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the 1st row
      {
         SMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 1UL, 0UL ) );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 4UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 8UL );

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

         if( submatrix(0,0) != 4 || submatrix(0,1) != 5 || submatrix(0,2) != -6 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 9 || submatrix(1,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 4  5 -6 )\n( 0  9  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the 1st row
      {
         SMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 1UL, 1UL ) );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 3UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 7UL );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) != 4 || submatrix(0,1) != 5 || submatrix(0,2) != -6 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 4  5 -6 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         SMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 1UL, 2UL ) );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 3UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 7UL );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) != 4 || submatrix(0,1) != 5 || submatrix(0,2) != -6 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 3UL, 0UL, 2UL, 4UL );

      // Erasing the 0th row
      {
         SMT::Iterator pos = submatrix.erase( 0UL, submatrix.begin( 0UL ), submatrix.end( 0UL ) );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 4UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 7UL );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) !=  0 || submatrix(0,2) != 0 || submatrix(0,3) !=  0 ||
             submatrix(1,0) != 7 || submatrix(1,1) != -8 || submatrix(1,2) != 9 || submatrix(1,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th row failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0  0  0  0 )\n( 7 -8  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 1st row
      {
         SMT::Iterator pos = submatrix.erase( 1UL, submatrix.begin( 1UL ), submatrix.find( 1UL, 2UL ) );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 2UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 5UL );

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

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 || submatrix(0,2) != 0 || submatrix(0,3) !=  0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) != 9 || submatrix(1,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the first half of the 1st row failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0  0  0  0 )\n( 0  0  9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 1st row
      {
         SMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 1UL, 2UL ), submatrix.end( 1UL ) );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 0UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 3UL );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 || submatrix(0,2) != 0 || submatrix(0,3) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) != 0 || submatrix(1,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the second half of the 1st row failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         SMT::Iterator pos = submatrix.erase( 1UL, submatrix.begin( 1UL ), submatrix.begin( 1UL ) );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 0UL );
         checkRows    ( mat_     , 5UL );
         checkColumns ( mat_     , 4UL );
         checkNonZeros( mat_     , 3UL );

         if( pos != submatrix.begin( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 || submatrix(0,2) != 0 || submatrix(0,3) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 || submatrix(1,2) != 0 || submatrix(1,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 1UL, 3UL, 3UL, 2UL );

      // Erasing the non-zero element at the end of the 1st column
      submatrix.erase( 2UL, 1UL );

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 5UL );
      checkRows    ( tmat_    , 4UL );
      checkColumns ( tmat_    , 5UL );
      checkNonZeros( tmat_    , 9UL );

      if( submatrix(0,0) !=  4 || submatrix(0,1) != -8 ||
          submatrix(1,0) !=  5 || submatrix(1,1) !=  9 ||
          submatrix(2,0) != -6 || submatrix(2,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n(  4 -8 )\n(  5  9 )\n( -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the 1st column
      submatrix.erase( 0UL, 1UL );

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 4UL );
      checkRows    ( tmat_    , 4UL );
      checkColumns ( tmat_    , 5UL );
      checkNonZeros( tmat_    , 8UL );

      if( submatrix(0,0) !=  4 || submatrix(0,1) != 0 ||
          submatrix(1,0) !=  5 || submatrix(1,1) != 9 ||
          submatrix(2,0) != -6 || submatrix(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n(  4 0 )\n(  5 9 )\n( -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the 1st column
      submatrix.erase( 1UL, 1UL );

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 3UL );
      checkRows    ( tmat_    , 4UL );
      checkColumns ( tmat_    , 5UL );
      checkNonZeros( tmat_    , 7UL );

      if( submatrix(0,0) !=  4 || submatrix(0,1) != 0 ||
          submatrix(1,0) !=  5 || submatrix(1,1) != 0 ||
          submatrix(2,0) != -6 || submatrix(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n(  4 0 )\n(  5 0 )\n( -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      submatrix.erase( 2UL, 1UL );

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 3UL );
      checkRows    ( tmat_    , 4UL );
      checkColumns ( tmat_    , 5UL );
      checkNonZeros( tmat_    , 7UL );

      if( submatrix(0,0) !=  4 || submatrix(0,1) != 0 ||
          submatrix(1,0) !=  5 || submatrix(1,1) != 0 ||
          submatrix(2,0) != -6 || submatrix(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 1UL, 3UL, 3UL, 2UL );

      // Erasing the non-zero element at the end of the 1st column
      {
         TSMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 2UL, 1UL ) );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 5UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 9UL );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) !=  4 || submatrix(0,1) != -8 ||
             submatrix(1,0) !=  5 || submatrix(1,1) !=  9 ||
             submatrix(2,0) != -6 || submatrix(2,1) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n(  4 -8 )\n(  5  9 )\n( -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the 1st column
      {
         TSMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 0UL, 1UL ) );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 4UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 8UL );

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

         if( submatrix(0,0) !=  4 || submatrix(0,1) != 0 ||
             submatrix(1,0) !=  5 || submatrix(1,1) != 9 ||
             submatrix(2,0) != -6 || submatrix(2,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n(  4 0 )\n(  5 9 )\n( -6 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the 1st column
      {
         TSMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 1UL, 1UL ) );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 3UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 7UL );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) !=  4 || submatrix(0,1) != 0 ||
             submatrix(1,0) !=  5 || submatrix(1,1) != 0 ||
             submatrix(2,0) != -6 || submatrix(2,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n(  4 0 )\n(  5 0 )\n( -6 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         TSMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 2UL, 1UL ) );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 3UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 7UL );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) !=  4 || submatrix(0,1) != 0 ||
             submatrix(1,0) !=  5 || submatrix(1,1) != 0 ||
             submatrix(2,0) != -6 || submatrix(2,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      TSMT submatrix = sub( tmat_, 0UL, 3UL, 4UL, 2UL );

      // Erasing the 0th column
      {
         TSMT::Iterator pos = submatrix.erase( 0UL, submatrix.begin( 0UL ), submatrix.end( 0UL ) );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 4UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 7UL );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) !=  7 ||
             submatrix(1,0) != 0 || submatrix(1,1) != -8 ||
             submatrix(2,0) != 0 || submatrix(2,1) !=  9 ||
             submatrix(3,0) != 0 || submatrix(3,1) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0  7 )\n( 0 -8 )\n( 0  9 )\n( 0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 1st column
      {
         TSMT::Iterator pos = submatrix.erase( 1UL, submatrix.begin( 1UL ), submatrix.find( 2UL, 1UL ) );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 2UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 5UL );

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

         if( submatrix(0,0) != 0 || submatrix(0,1) !=  0 ||
             submatrix(1,0) != 0 || submatrix(1,1) !=  0 ||
             submatrix(2,0) != 0 || submatrix(2,1) !=  9 ||
             submatrix(3,0) != 0 || submatrix(3,1) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0  0 )\n( 0  0 )\n( 0  9 )\n( 0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 1st column
      {
         TSMT::Iterator pos = submatrix.erase( 1UL, submatrix.find( 2UL, 1UL ), submatrix.end( 1UL ) );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 3UL );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 ||
             submatrix(2,0) != 0 || submatrix(2,1) != 0 ||
             submatrix(3,0) != 0 || submatrix(3,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
                << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         TSMT::Iterator pos = submatrix.erase( 1UL, submatrix.begin( 1UL ), submatrix.begin( 1UL ) );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkRows    ( tmat_    , 4UL );
         checkColumns ( tmat_    , 5UL );
         checkNonZeros( tmat_    , 3UL );

         if( pos != submatrix.begin( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( submatrix(0,0) != 0 || submatrix(0,1) != 0 ||
             submatrix(1,0) != 0 || submatrix(1,1) != 0 ||
             submatrix(2,0) != 0 || submatrix(2,1) != 0 ||
             submatrix(3,0) != 0 || submatrix(3,1) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the 0th column failed\n"
                << " Details:\n"
                << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat, 1UL, 0UL, 1UL, 20UL );

      // Increasing the capacity of the row
      submatrix.reserve( 10UL );

      checkRows    ( submatrix,  1UL );
      checkColumns ( submatrix, 20UL );
      checkCapacity( submatrix, 10UL );
      checkNonZeros( submatrix,  0UL );

      // Further increasing the capacity of the row
      submatrix.reserve( 15UL );

      checkRows    ( submatrix,  1UL );
      checkColumns ( submatrix, 20UL );
      checkCapacity( submatrix, 15UL );
      checkNonZeros( submatrix,  0UL );
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Columnt-major SparseSubmatrix::reserve()";

      TMT mat( 20UL, 3UL );

      TSMT submatrix = sub( mat, 0UL, 1UL, 20UL, 1UL );

      // Increasing the capacity of the column
      submatrix.reserve( 10UL );

      checkRows    ( submatrix, 20UL );
      checkColumns ( submatrix,  1UL );
      checkCapacity( submatrix, 10UL );
      checkNonZeros( submatrix,  0UL );

      // Further increasing the capacity of the column
      submatrix.reserve( 15UL );

      checkRows    ( submatrix, 20UL );
      checkColumns ( submatrix,  1UL );
      checkCapacity( submatrix, 15UL );
      checkNonZeros( submatrix,  0UL );
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
      SMT submatrix = sub( mat_, 2UL, 1UL, 2UL, 2UL );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 3UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 2UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != -3 ||
          submatrix(1,0) != 4 || submatrix(1,1) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 0 -3 )\n( 4  5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      submatrix.scale( 2 );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 3UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 2UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != -6 ||
          submatrix(1,0) != 8 || submatrix(1,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n( 0 -6 )\n( 8 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      submatrix.scale( 0.5 );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 3UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 2UL );

      if( submatrix(0,0) != 0 || submatrix(0,1) != -3 ||
          submatrix(1,0) != 4 || submatrix(1,1) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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
      TSMT submatrix = sub( tmat_, 1UL, 2UL, 2UL, 2UL );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 3UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 2UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 4 ||
          submatrix(1,0) != -3 || submatrix(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n(  0 4 )\n( -3 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      submatrix.scale( 2 );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 3UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 2UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) !=  8 ||
          submatrix(1,0) != -6 || submatrix(1,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
             << "   Expected result:\n(  0  8 )\n( -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      submatrix.scale( 0.5 );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 3UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 2UL );

      if( submatrix(0,0) !=  0 || submatrix(0,1) != 4 ||
          submatrix(1,0) != -3 || submatrix(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << submatrix << "\n"
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

      SMT submatrix = sub( mat_, 1UL, 1UL, 3UL, 2UL );

      checkRows    ( submatrix, 3UL );
      checkColumns ( submatrix, 2UL );
      checkNonZeros( submatrix, 4UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 1UL );
      checkNonZeros( submatrix, 2UL, 2UL );

      // Searching for the first element
      {
         ConstIterator pos( submatrix.find( 0UL, 0UL ) );

         if( pos == submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( submatrix.find( 1UL, 1UL ) );

         if( pos == submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( submatrix.find( 1UL, 0UL ) );

         if( pos != submatrix.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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

      TSMT submatrix = sub( tmat_, 1UL, 1UL, 2UL, 3UL );

      checkRows    ( submatrix, 2UL );
      checkColumns ( submatrix, 3UL );
      checkNonZeros( submatrix, 4UL );
      checkNonZeros( submatrix, 0UL, 1UL );
      checkNonZeros( submatrix, 1UL, 1UL );
      checkNonZeros( submatrix, 2UL, 2UL );

      // Searching for the first element
      {
         ConstIterator pos( submatrix.find( 0UL, 0UL ) );

         if( pos == submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( submatrix.find( 1UL, 2UL ) );

         if( pos == submatrix.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( submatrix.find( 1UL, 0UL ) );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 1UL, 4UL );

      checkRows    ( submatrix, 1UL );
      checkColumns ( submatrix, 4UL );
      checkNonZeros( submatrix, 1UL );
      checkNonZeros( submatrix, 0UL, 1UL );

      // Determining the lower bound for position (0,0)
      {
         ConstIterator pos( submatrix.lowerBound( 0UL, 0UL ) );

         if( pos == submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (0,1)
      {
         ConstIterator pos( submatrix.lowerBound( 0UL, 1UL ) );

         if( pos == submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (0,2)
      {
         ConstIterator pos( submatrix.lowerBound( 0UL, 2UL ) );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,2)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 4UL, 1UL );

      checkRows    ( submatrix, 4UL );
      checkColumns ( submatrix, 1UL );
      checkNonZeros( submatrix, 1UL );
      checkNonZeros( submatrix, 0UL, 1UL );

      // Determining the lower bound for position (0,0)
      {
         ConstIterator pos( submatrix.lowerBound( 0UL, 0UL ) );

         if( pos == submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,0)
      {
         ConstIterator pos( submatrix.lowerBound( 1UL, 0UL ) );

         if( pos == submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,0)
      {
         ConstIterator pos( submatrix.lowerBound( 2UL, 0UL ) );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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

      SMT submatrix = sub( mat_, 1UL, 0UL, 1UL, 4UL );

      checkRows    ( submatrix, 1UL );
      checkColumns ( submatrix, 4UL );
      checkNonZeros( submatrix, 1UL );
      checkNonZeros( submatrix, 0UL, 1UL );

      // Determining the upper bound for position (0,0)
      {
         ConstIterator pos( submatrix.upperBound( 0UL, 0UL ) );

         if( pos == submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (0,1)
      {
         ConstIterator pos( submatrix.upperBound( 0UL, 1UL ) );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (0,2)
      {
         ConstIterator pos( submatrix.upperBound( 0UL, 2UL ) );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,2)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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

      TSMT submatrix = sub( tmat_, 0UL, 1UL, 4UL, 1UL );

      checkRows    ( submatrix, 4UL );
      checkColumns ( submatrix, 1UL );
      checkNonZeros( submatrix, 1UL );
      checkNonZeros( submatrix, 0UL, 1UL );

      // Determining the upper bound for position (0,0)
      {
         ConstIterator pos( submatrix.upperBound( 0UL, 0UL ) );

         if( pos == submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,0)
      {
         ConstIterator pos( submatrix.upperBound( 1UL, 0UL ) );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,0)
      {
         ConstIterator pos( submatrix.upperBound( 2UL, 0UL ) );

         if( pos != submatrix.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,0)\n"
                << "   Current submatrix:\n" << submatrix << "\n";
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
         SMT submatrix = sub( mat_, 0UL, 0UL, 1UL, 4UL );

         if( isDefault( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default submatrix
      {
         SMT submatrix = sub( mat_, 1UL, 0UL, 1UL, 4UL );

         if( isDefault( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
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
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 4UL, 1UL );

         if( isDefault( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default submatrix
      {
         TSMT submatrix = sub( tmat_, 0UL, 1UL, 4UL, 1UL );

         if( isDefault( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
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
      test_ = "Row-major isnan()";

      typedef blaze::CompressedMatrix<float,blaze::rowMajor>  MatrixType;
      typedef blaze::SparseSubmatrix<MatrixType>              SubmatrixType;

      MatrixType mat( mat_ );

      // isnan with empty 2x2 matrix
      {
         SubmatrixType submatrix = sub( mat, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );

         if( isnan( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with filled 2x3 matrix
      {
         SubmatrixType submatrix = sub( mat, 2UL, 1UL, 2UL, 3UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 4UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 3UL );

         if( isnan( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major isnan()";

      typedef blaze::CompressedMatrix<float,blaze::columnMajor>  MatrixType;
      typedef blaze::SparseSubmatrix<MatrixType>                 SubmatrixType;

      MatrixType mat( tmat_ );

      // isnan with empty 2x2 matrix
      {
         SubmatrixType submatrix = sub( mat, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );

         if( isnan( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with filled 3x2 matrix
      {
         SubmatrixType submatrix = sub( mat, 1UL, 2UL, 3UL, 2UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 4UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 3UL );

         if( isnan( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
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
      test_ = "Row-major isDiagonal()";

      initialize();
      mat_(0,0) = 11;
      mat_(2,0) =  0;

      // Non-quadratic submatrix
      {
         SMT submatrix = sub( mat_, 0UL, 0UL, 2UL, 3UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );

         if( isDiagonal( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         SMT submatrix = sub( mat_, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );

         if( isDiagonal( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         SMT submatrix = sub( mat_, 0UL, 0UL, 3UL, 3UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );
         checkNonZeros( submatrix, 2UL, 1UL );

         if( isDiagonal( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-diagonal matrix
      {
         SMT submatrix = sub( mat_, 0UL, 0UL, 4UL, 4UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 6UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );
         checkNonZeros( submatrix, 2UL, 1UL );
         checkNonZeros( submatrix, 3UL, 3UL );

         if( isDiagonal( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDiagonal()";

      initialize();
      tmat_(0,0) = 11;
      tmat_(0,2) =  0;

      // Non-quadratic submatrix
      {
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 3UL, 2UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );

         if( isDiagonal( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         TSMT submatrix = sub( tmat_, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );

         if( isDiagonal( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 3UL, 3UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );
         checkNonZeros( submatrix, 2UL, 1UL );

         if( isDiagonal( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-diagonal matrix
      {
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 4UL, 4UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 6UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );
         checkNonZeros( submatrix, 2UL, 1UL );
         checkNonZeros( submatrix, 3UL, 3UL );

         if( isDiagonal( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
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
      test_ = "Row-major isSymmetric()";

      initialize();
      mat_(0,0) = 11;
      mat_(2,0) =  0;
      mat_(2,3) =  5;
      mat_(3,1) =  0;

      // Non-quadratic matrix
      {
         SMT submatrix = sub( mat_, 0UL, 0UL, 2UL, 3UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );

         if( isSymmetric( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         SMT submatrix = sub( mat_, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );

         if( isSymmetric( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         SMT submatrix = sub( mat_, 0UL, 0UL, 3UL, 3UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );
         checkNonZeros( submatrix, 2UL, 1UL );

         if( isSymmetric( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix
      {
         SMT submatrix = sub( mat_, 1UL, 0UL, 4UL, 4UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 9UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 2UL );
         checkNonZeros( submatrix, 2UL, 2UL );
         checkNonZeros( submatrix, 3UL, 4UL );

         if( isSymmetric( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Symmetric matrix
      {
         SMT submatrix = sub( mat_, 0UL, 0UL, 4UL, 4UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 6UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );
         checkNonZeros( submatrix, 2UL, 2UL );
         checkNonZeros( submatrix, 3UL, 2UL );

         if( isSymmetric( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major isSymmetric()";

      initialize();
      tmat_(0,0) = 11;
      tmat_(0,2) =  0;
      tmat_(3,2) =  5;
      tmat_(1,3) =  0;

      // Non-quadratic matrix
      {
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 3UL, 2UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );

         if( isSymmetric( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         TSMT submatrix = sub( tmat_, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );
         checkNonZeros( submatrix, 0UL, 0UL );
         checkNonZeros( submatrix, 1UL, 0UL );

         if( isSymmetric( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 3UL, 3UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 3UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );
         checkNonZeros( submatrix, 2UL, 1UL );

         if( isSymmetric( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix
      {
         TSMT submatrix = sub( tmat_, 0UL, 1UL, 4UL, 4UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 9UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 2UL );
         checkNonZeros( submatrix, 2UL, 2UL );
         checkNonZeros( submatrix, 3UL, 4UL );

         if( isSymmetric( submatrix ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Symmetric matrix
      {
         TSMT submatrix = sub( tmat_, 0UL, 0UL, 4UL, 4UL );

         checkRows    ( submatrix, 4UL );
         checkColumns ( submatrix, 4UL );
         checkNonZeros( submatrix, 6UL );
         checkNonZeros( submatrix, 0UL, 1UL );
         checkNonZeros( submatrix, 1UL, 1UL );
         checkNonZeros( submatrix, 2UL, 2UL );
         checkNonZeros( submatrix, 3UL, 2UL );

         if( isSymmetric( submatrix ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Submatrix:\n" << submatrix << "\n";
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
      test_ = "Row-major min()";

      initialize();

      // Attempt to find the minimum in an empty submatrix
      {
         SMT submatrix = sub( mat_, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );

         const int minimum = min( submatrix );

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
         SMT submatrix = sub( mat_, 1UL, 1UL, 2UL, 3UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 2UL );

         const int minimum = min( submatrix );

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
         SMT submatrix = sub( mat_, 3UL, 1UL, 2UL, 3UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 6UL );

         const int minimum = min( submatrix );

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
      test_ = "Column-major min()";

      initialize();

      // Attempt to find the minimum in an empty submatrix
      {
         TSMT submatrix = sub( tmat_, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );

         const int minimum = min( submatrix );

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
         TSMT submatrix = sub( tmat_, 1UL, 1UL, 3UL, 2UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 2UL );

         const int minimum = min( submatrix );

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
         TSMT submatrix = sub( tmat_, 1UL, 3UL, 3UL, 2UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 6UL );

         const int minimum = min( submatrix );

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
      test_ = "Row-major max()";

      initialize();

      // Attempt to find the maximum in an empty submatrix
      {
         SMT submatrix = sub( mat_, 0UL, 2UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );

         const int maximum = max( submatrix );

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
         SMT submatrix = sub( mat_, 1UL, 1UL, 2UL, 3UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 2UL );

         const int maximum = max( submatrix );

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
         SMT submatrix = sub( mat_, 3UL, 1UL, 2UL, 3UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 3UL );
         checkNonZeros( submatrix, 6UL );

         const int maximum = max( submatrix );

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
      test_ = "Column-major max()";

      initialize();

      // Attempt to find the maximum in an empty submatrix
      {
         TSMT submatrix = sub( tmat_, 2UL, 0UL, 2UL, 2UL );

         checkRows    ( submatrix, 2UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 0UL );

         const int maximum = max( submatrix );

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
         TSMT submatrix = sub( tmat_, 1UL, 1UL, 3UL, 2UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 2UL );

         const int maximum = max( submatrix );

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
         TSMT submatrix = sub( tmat_, 1UL, 3UL, 3UL, 2UL );

         checkRows    ( submatrix, 3UL );
         checkColumns ( submatrix, 2UL );
         checkNonZeros( submatrix, 6UL );

         const int maximum = max( submatrix );

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
