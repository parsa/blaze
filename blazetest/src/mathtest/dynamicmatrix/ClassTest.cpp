//=================================================================================================
/*!
//  \file src/mathtest/dynamicmatrix/ClassTest.cpp
//  \brief Source file for the DynamicMatrix class test
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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/dynamicmatrix/ClassTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace mathtest {

namespace dynamicmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DynamicMatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testAlignment< signed char          >( "signed char"          );
   testAlignment< unsigned char        >( "unsigned char"        );
   testAlignment< short                >( "short"                );
   testAlignment< unsigned short       >( "unsigned short"       );
   testAlignment< int                  >( "int"                  );
   testAlignment< unsigned int         >( "unsigned int"         );
   testAlignment< float                >( "float"                );
   testAlignment< double               >( "double"               );
   testAlignment< long double          >( "long double"          );
   testAlignment< complex<float>       >( "complex<float>"       );
   testAlignment< complex<double>      >( "complex<double>"      );
   testAlignment< complex<long double> >( "complex<long double>" );

   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testDivAssign();
   testFunctionCall();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testExtend();
   testReserve();
   testTranspose();
   testIsDiagonal();
   testIsSymmetric();
   testScale();
   testSwap();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the DynamicMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the DynamicMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor
   {
      test_ = "Row-major DynamicMatrix default constructor";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix size constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix size constructor (0x4)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 4UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix size constructor (3x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix size constructor (3x4)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
   }


   //=====================================================================================
   // Row-major homogeneous initialization
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix homogeneous initialization constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 0UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix homogeneous initialization constructor (0x4)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 4UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix homogeneous initialization constructor (3x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL, 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix homogeneous initialization constructor (3x4)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 4UL );
      checkNonZeros( mat,  1UL, 4UL );
      checkNonZeros( mat,  2UL, 4UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array initialization
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix array initialization constructor";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::rowMajor> mat( array );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix copy constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 0UL, 0UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix copy constructor (0x3)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 0UL, 3UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix copy constructor (2x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 0UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix copy constructor (2x3)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

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
      test_ = "Column-major DynamicMatrix default constructor";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix size constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix size constructor (0x4)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 4UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix size constructor (3x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix size constructor (3x4)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 4UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
   }


   //=====================================================================================
   // Column-major homogeneous initialization
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix homogeneous initialization constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 0UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix homogeneous initialization constructor (0x4)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 4UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix homogeneous initialization constructor (3x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL, 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix homogeneous initialization constructor (3x4)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 4UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array initialization
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix array initialization constructor";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::columnMajor> mat( array );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix copy constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 0UL, 0UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix copy constructor (0x3)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 0UL, 3UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix copy constructor (2x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 0UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix copy constructor (2x3)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
/*!\brief Test of the DynamicMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the DynamicMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix homogeneous assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
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
      test_ = "Row-major DynamicMatrix array assignment";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      mat = array;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

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
      test_ = "Row-major DynamicMatrix copy assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

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
      test_ = "Row-major DynamicMatrix copy assignment stress test";

      typedef blaze::DynamicMatrix<int,blaze::rowMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment";

      blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat1;
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

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
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment stress test";

      typedef blaze::DynamicMatrix<unsigned int,blaze::rowMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const unsigned int min( randmin );
      const unsigned int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

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
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment stress test";

      typedef blaze::DynamicMatrix<unsigned int,blaze::columnMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const unsigned int min( randmin );
      const unsigned int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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
      test_ = "Row-major/row-major DynamicMatrix sparse matrix assignment stress test";

      typedef blaze::CompressedMatrix<int,blaze::rowMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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
      test_ = "Row-major/column-major DynamicMatrix sparse matrix assignment stress test";

      typedef blaze::CompressedMatrix<int,blaze::columnMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix homogeneous assigment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 4UL );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
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
      test_ = "Column-major DynamicMatrix array initialization constructor";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      mat = array;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

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
      test_ = "Column-major DynamicMatrix copy assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
      test_ = "Row-major DynamicMatrix copy assignment stress test";

      typedef blaze::DynamicMatrix<int,blaze::columnMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment stress test";

      typedef blaze::DynamicMatrix<unsigned int,blaze::rowMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const unsigned int min( randmin );
      const unsigned int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment";

      blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat1;
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment stress test";

      typedef blaze::DynamicMatrix<unsigned int,blaze::columnMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const unsigned int min( randmin );
      const unsigned int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

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
      test_ = "Column-major/row-major DynamicMatrix sparse matrix assignment stress test";

      typedef blaze::CompressedMatrix<int,blaze::rowMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

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
      test_ = "Column-major/column-major DynamicMatrix sparse matrix assignment stress test";

      typedef blaze::CompressedMatrix<int,blaze::columnMajor>  RandomMatrixType;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

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
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the DynamicMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix addition assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 0 );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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
      test_ = "Row-major/column-major DynamicMatrix dense matrix addition assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 0 );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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
      test_ = "Row-major/column-major DynamicMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix addition assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 0 );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
      test_ = "Column-major/column-major DynamicMatrix dense matrix addition assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 0 );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

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


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
      test_ = "Column-major/column-major DynamicMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the DynamicMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix subtraction assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 0 );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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
      test_ = "Row-major/column-major DynamicMatrix dense matrix subtraction assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 0 );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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
      test_ = "Row-major/column-major DynamicMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

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


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix subtraction assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 0 );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
      test_ = "Column-major/column-major DynamicMatrix dense matrix subtraction assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 0 );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

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


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
      test_ = "Column-major/column-major DynamicMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

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
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the DynamicMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Row-major scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major scalar multiplication assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
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
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major scalar multiplication assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
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
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicMatrix division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the DynamicMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testDivAssign()
{
   //=====================================================================================
   // Row-major scalar division assignment
   //=====================================================================================

   {
      test_ = "Row-major scalar division assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
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
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar division assignment
   //=====================================================================================

   {
      test_ = "Column-major scalar division assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
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
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the DynamicMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::operator()";

      // Writing the first element
      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 0 );
      mat(2,1) = 1;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
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

      // Writing the second element
      mat(1,4) = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat,  2UL );
      checkNonZeros( mat,  0UL, 0UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat(2,1) != 1 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the third element
      mat(0,3) = 3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat,  3UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 1UL );

      if( mat(2,1) != 1 || mat(1,4) != 2 || mat(0,3) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the fourth element
      mat(2,2) = 4;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat,  4UL );
      checkNonZeros( mat,  0UL, 1UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(2,1) != 1 || mat(1,4) != 2 || mat(0,3) != 3 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::operator()";

      // Writing the first element
      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 0 );
      mat(2,1) = 1;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
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

      // Writing the second element
      mat(1,4) = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
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

      // Writing the third element
      mat(0,3) = 3;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
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

      // Writing the fourth element
      mat(2,2) = 4;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
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
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the nonZeros member function of DynamicMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of DynamicMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::nonZeros()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

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
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
         mat(0,1) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 0;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );

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
      test_ = "Column-major DynamicMatrix::nonZeros()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

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
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );
         mat(0,1) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 0;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

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
/*!\brief Test of the reset member function of DynamicMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of DynamicMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::reset()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(0,2) = 3;
      mat(1,0) = 4;
      mat(1,1) = 5;
      mat(1,2) = 6;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

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

      // Resetting row 1
      mat.reset( 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 0UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      mat.reset();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );

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
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::reset()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(0,2) = 3;
      mat(1,0) = 4;
      mat(1,1) = 5;
      mat(1,2) = 6;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

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

      // Resetting column 1
      mat.reset( 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      mat.reset();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );

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
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the clear member function of DynamicMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the clear member function of DynamicMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::clear()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(0,2) = 3;
      mat(1,0) = 4;
      mat(1,1) = 5;
      mat(1,2) = 6;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

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

      // Clearing the matrix
      mat.clear();

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::clear()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(0,2) = 3;
      mat(1,0) = 4;
      mat(1,1) = 5;
      mat(1,2) = 6;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

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

      // Clearing the matrix
      mat.clear();

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the resize member function of DynamicMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the resize member function of DynamicMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::resize()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 0x3
      mat.resize( 0UL, 3UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 5x0
      mat.resize( 5UL, 0UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 2x1
      mat.resize( 2UL, 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 2UL );

      // Resizing to 3x2 and preserving the elements
      mat(0,0) = 1;
      mat(1,0) = 2;
      mat.resize( 3UL, 2UL, true );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 x )\n( 2 x )\n( x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2 and preserving the elements
      mat(0,1) = 3;
      mat(1,1) = 4;
      mat.resize( 2UL, 2UL, true );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 3 || mat(1,0) != 2 || mat(1,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 3 )\n( 2 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 1x1
      mat.resize( 1UL, 1UL );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 1UL );

      // Resizing to 0x0
      mat.resize( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::resize()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 0x3
      mat.resize( 0UL, 3UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 5x0
      mat.resize( 5UL, 0UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 2x1
      mat.resize( 2UL, 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 2UL );

      // Resizing to 3x2 and preserving the elements
      mat(0,0) = 1;
      mat(1,0) = 2;
      mat.resize( 3UL, 2UL, true );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 x )\n( 2 x )\n( x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2 and preserving the elements
      mat(0,1) = 3;
      mat(1,1) = 4;
      mat.resize( 2UL, 2UL, true );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 3 || mat(1,0) != 2 || mat(1,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 3 )\n( 2 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 1x1
      mat.resize( 1UL, 1UL );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 1UL );

      // Resizing to 0x0
      mat.resize( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the extend member function of DynamicMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the extend member function of DynamicMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::extend()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Increasing the size of the matrix
      mat.extend( 2UL, 2UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );

      // Further increasing the size of the matrix and preserving the elements
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(1,0) = 3;
      mat(1,1) = 4;
      mat.extend( 1UL, 1UL, true );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 9UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 ||
          mat(1,0) != 3 || mat(1,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 x )\n( 3 4 x )\n( x x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Further increasing the size of the matrix
      mat.extend( 4UL, 10UL, false );

      checkRows    ( mat,  7UL );
      checkColumns ( mat, 13UL );
      checkCapacity( mat, 91UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::extend()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Increasing the size of the matrix
      mat.extend( 2UL, 2UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );

      // Further increasing the size of the matrix and preserving the elements
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(1,0) = 3;
      mat(1,1) = 4;
      mat.extend( 1UL, 1UL, true );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 9UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 ||
          mat(1,0) != 3 || mat(1,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 x )\n( 3 4 x )\n( x x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Further increasing the size of the matrix
      mat.extend( 4UL, 10UL, false );

      checkRows    ( mat,  7UL );
      checkColumns ( mat, 13UL );
      checkCapacity( mat, 91UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reserve member function of DynamicMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reserve member function of DynamicMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::reserve()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Increasing the capacity of the matrix
      mat.reserve( 10UL );

      checkRows    ( mat,  0UL );
      checkColumns ( mat,  0UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat,  0UL );

      // Further increasing the capacity of the matrix
      mat.reserve( 20UL );

      checkRows    ( mat,  0UL );
      checkColumns ( mat,  0UL );
      checkCapacity( mat, 20UL );
      checkNonZeros( mat,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::reserve()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Increasing the capacity of the matrix
      mat.reserve( 10UL );

      checkRows    ( mat,  0UL );
      checkColumns ( mat,  0UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat,  0UL );

      // Further increasing the capacity of the matrix
      mat.reserve( 20UL );

      checkRows    ( mat,  0UL );
      checkColumns ( mat,  0UL );
      checkCapacity( mat, 20UL );
      checkNonZeros( mat,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the transpose member function of the DynamicMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the transpose member functions of the DynamicMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::transpose()";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 0 );
      mat(0,0) = 1;
      mat(0,2) = 2;
      mat(0,4) = 3;
      mat(1,1) = 4;
      mat(1,3) = 5;
      mat(2,0) = 6;
      mat(2,2) = 7;
      mat(2,4) = 8;

      mat.transpose();

      checkRows    ( mat,  5UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat,  8UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 1UL );
      checkNonZeros( mat,  2UL, 2UL );
      checkNonZeros( mat,  3UL, 1UL );
      checkNonZeros( mat,  4UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 6 ||
          mat(1,0) != 0 || mat(1,1) != 4 || mat(1,2) != 0 ||
          mat(2,0) != 2 || mat(2,1) != 0 || mat(2,2) != 7 ||
          mat(3,0) != 0 || mat(3,1) != 5 || mat(3,2) != 0 ||
          mat(4,0) != 3 || mat(4,1) != 0 || mat(4,2) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 0 4 0 )\n( 2 0 7 )\n( 0 5 0 )\n( 3 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::transpose()";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 0 );
      mat(0,0) = 1;
      mat(0,2) = 2;
      mat(0,4) = 3;
      mat(1,1) = 4;
      mat(1,3) = 5;
      mat(2,0) = 6;
      mat(2,2) = 7;
      mat(2,4) = 8;

      mat.transpose();

      checkRows    ( mat,  5UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat,  8UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 6 ||
          mat(1,0) != 0 || mat(1,1) != 4 || mat(1,2) != 0 ||
          mat(2,0) != 2 || mat(2,1) != 0 || mat(2,2) != 7 ||
          mat(3,0) != 0 || mat(3,1) != 5 || mat(3,2) != 0 ||
          mat(4,0) != 3 || mat(4,1) != 0 || mat(4,2) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 0 4 0 )\n( 2 0 7 )\n( 0 5 0 )\n( 3 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDiagonal member function of the DynamicMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDiagonal member functions of the DynamicMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDiagonal()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::isDiagonal()";

      // Non-quadratic matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( mat.isDiagonal() != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat.isDiagonal() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat.isDiagonal() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-diagonal matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat.isDiagonal() != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
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
      test_ = "Column-major DynamicMatrix::isDiagonal()";

      // Non-quadratic matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat.isDiagonal() != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat.isDiagonal() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat.isDiagonal() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-diagonal matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat.isDiagonal() != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isSymmetric member function of the DynamicMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isSymmetric member functions of the DynamicMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsSymmetric()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::isSymmetric()";

      // Non-quadratic matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( mat.isSymmetric() != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat.isSymmetric() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat.isSymmetric() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat.isSymmetric() != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Symmetric matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat.isSymmetric() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
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
      test_ = "Column-major DynamicMatrix::isSymmetric()";

      // Non-quadratic matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat.isSymmetric() != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat.isSymmetric() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat.isSymmetric() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat.isSymmetric() != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Symmetric matrix
      {
         // Initialization check
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat.isSymmetric() != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of DynamicMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of DynamicMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScale()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::scale()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL );
      mat(0,0) = 1;
      mat(0,1) = 2;
      mat(1,0) = 3;
      mat(1,1) = 4;
      mat(2,0) = 5;
      mat(2,1) = 6;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

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

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

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

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

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
      using blaze::complex;

      blaze::DynamicMatrix<complex<float>,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) = complex<float>( 1.0F, 0.0F );
      mat(0,1) = complex<float>( 2.0F, 0.0F );
      mat(1,0) = complex<float>( 3.0F, 0.0F );
      mat(1,1) = complex<float>( 4.0F, 0.0F );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );

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
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::scale()";

      // Initialization check
      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL );
      mat(0,0) = 1;
      mat(0,1) = 4;
      mat(1,0) = 2;
      mat(1,1) = 5;
      mat(2,0) = 3;
      mat(2,1) = 6;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

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

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

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

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

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
      using blaze::complex;

      blaze::DynamicMatrix<complex<float>,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) = complex<float>( 1.0F, 0.0F );
      mat(0,1) = complex<float>( 2.0F, 0.0F );
      mat(1,0) = complex<float>( 3.0F, 0.0F );
      mat(1,1) = complex<float>( 4.0F, 0.0F );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );

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
/*!\brief Test of the swap functionality of the DynamicMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the swap function of the DynamicMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix swap";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 2UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 0;
      mat1(1,1) = 3;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 2UL );
      mat2(0,0) = 4;
      mat2(0,1) = 3;
      mat2(1,0) = 2;
      mat2(1,1) = 1;

      swap( mat1, mat2 );

      checkRows    ( mat1, 2UL );
      checkColumns ( mat1, 2UL );
      checkCapacity( mat1, 4UL );
      checkNonZeros( mat1, 4UL );
      checkNonZeros( mat1, 0UL, 2UL );
      checkNonZeros( mat1, 1UL, 2UL );

      if( mat1(0,0) != 4 || mat1(0,1) != 3 || mat1(1,0) != 2 || mat1(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n( 4 3 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 2UL );
      checkCapacity( mat2, 4UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );

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
      test_ = "Column-major DynamicMatrix swap";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 2UL );
      mat1(0,0) = 1;
      mat1(0,1) = 0;
      mat1(1,0) = 2;
      mat1(1,1) = 3;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 2UL, 2UL );
      mat2(0,0) = 4;
      mat2(0,1) = 2;
      mat2(1,0) = 3;
      mat2(1,1) = 1;

      swap( mat1, mat2 );

      checkRows    ( mat1, 2UL );
      checkColumns ( mat1, 2UL );
      checkCapacity( mat1, 4UL );
      checkNonZeros( mat1, 4UL );
      checkNonZeros( mat1, 0UL, 2UL );
      checkNonZeros( mat1, 1UL, 2UL );

      if( mat1(0,0) != 4 || mat1(0,1) != 2 || mat1(1,0) != 3 || mat1(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n( 4 2 )\n( 3 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 2UL );
      checkCapacity( mat2, 4UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );

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

} // namespace dynamicmatrix

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
   std::cout << "   Running DynamicMatrix class test..." << std::endl;

   try
   {
      RUN_DYNAMICMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during DynamicMatrix class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
