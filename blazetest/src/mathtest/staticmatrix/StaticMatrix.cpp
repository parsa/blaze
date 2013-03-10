//=================================================================================================
/*!
//  \file src/mathtest/staticmatrix/StaticMatrix.cpp
//  \brief Source file for the StaticMatrix test
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
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/StaticMatrix.h>


namespace blazetest {

namespace mathtest {

namespace staticmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the StaticMatrix test class.
//
// \exception std::runtime_error Operation error detected.
*/
StaticMatrix::StaticMatrix()
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
   testFunctionCall();
   testNonZeros();
   testReset();
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
/*!\brief Test of the StaticMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the StaticMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticMatrix::testConstructors()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   // Default constructor
   {
      test_ = "Row-major StaticMatrix default constructor";

      blaze::StaticMatrix<int,3UL,4UL,blaze::rowMajor> mat;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat,  0UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 || mat(1,3) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   // Homogeneous initialization
   {
      test_ = "Row-major StaticMatrix homogeneous initialization constructor";

      blaze::StaticMatrix<int,3UL,4UL,blaze::rowMajor> mat( 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );

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

   // 1x2 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x2 initialization constructor";

      blaze::StaticMatrix<int,1UL,2UL,blaze::rowMajor> mat( 1, 2 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 2x1 initialization constructor";

      blaze::StaticMatrix<int,2UL,1UL,blaze::rowMajor> mat( 1, 2 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x3 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x3 initialization constructor";

      blaze::StaticMatrix<int,1UL,3UL,blaze::rowMajor> mat( 1, 2, 3 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 3x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 3x1 initialization constructor";

      blaze::StaticMatrix<int,3UL,1UL,blaze::rowMajor> mat( 1, 2, 3 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x4 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x4 initialization constructor";

      blaze::StaticMatrix<int,1UL,4UL,blaze::rowMajor> mat( 1, 2, 3, 4 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x2 initialization constructor
   {
      test_ = "Row-major StaticMatrix 2x2 initialization constructor";

      blaze::StaticMatrix<int,2UL,2UL,blaze::rowMajor> mat( 1, 2, 3, 4 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(1,0) != 3 || mat(1,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 4x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 4x1 initialization constructor";

      blaze::StaticMatrix<int,4UL,1UL,blaze::rowMajor> mat( 1, 2, 3, 4 );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x5 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x5 initialization constructor";

      blaze::StaticMatrix<int,1UL,5UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 5UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 || mat(0,4) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 5x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 5x1 initialization constructor";

      blaze::StaticMatrix<int,5UL,1UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5 );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 5UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 || mat(4,0) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x6 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x6 initialization constructor";

      blaze::StaticMatrix<int,1UL,6UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 6UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(0,3) != 4 || mat(0,4) != 5 || mat(0,5) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x3 initialization constructor
   {
      test_ = "Row-major StaticMatrix 2x3 initialization constructor";

      blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

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

   // 3x2 initialization constructor
   {
      test_ = "Row-major StaticMatrix 3x2 initialization constructor";

      blaze::StaticMatrix<int,3UL,2UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 ||
          mat(1,0) != 3 || mat(1,1) != 4 ||
          mat(2,0) != 5 || mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 )\n( 3 4 )\n( 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 6x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 6x1 initialization constructor";

      blaze::StaticMatrix<int,6UL,1UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 ||
          mat(3,0) != 4 || mat(4,0) != 5 || mat(5,0) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n( 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x7 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x7 initialization constructor";

      blaze::StaticMatrix<int,1UL,7UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 7UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 ||
          mat(0,4) != 5 || mat(0,5) != 6 || mat(0,6) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 7x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 7x1 initialization constructor";

      blaze::StaticMatrix<int,7UL,1UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7 );

      checkRows    ( mat, 7UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 ||
          mat(4,0) != 5 || mat(5,0) != 6 || mat(6,0) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n( 6 )\n( 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x8 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x8 initialization constructor";

      blaze::StaticMatrix<int,1UL,8UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 8UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 8UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 ||
          mat(0,4) != 5 || mat(0,5) != 6 || mat(0,6) != 7 || mat(0,7) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x4 initialization constructor
   {
      test_ = "Row-major StaticMatrix 2x4 initialization constructor";

      blaze::StaticMatrix<int,2UL,4UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 8UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 ||
          mat(1,0) != 5 || mat(1,1) != 6 || mat(1,2) != 7 || mat(1,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n( 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 4x2 initialization constructor
   {
      test_ = "Row-major StaticMatrix 4x2 initialization constructor";

      blaze::StaticMatrix<int,4UL,2UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8 );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 8UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(1,0) != 3 || mat(1,1) != 4 ||
          mat(2,0) != 5 || mat(2,1) != 6 || mat(3,0) != 7 || mat(3,1) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 )\n( 3 4 )\n( 5 6 )\n( 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 8x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 8x1 initialization constructor";

      blaze::StaticMatrix<int,8UL,1UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8 );

      checkRows    ( mat, 8UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 8UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 ||
          mat(4,0) != 5 || mat(5,0) != 6 || mat(6,0) != 7 || mat(7,0) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n( 6 )\n( 7 )\n( 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x9 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x9 initialization constructor";

      blaze::StaticMatrix<int,1UL,9UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 9UL );
      checkCapacity( mat, 9UL );
      checkNonZeros( mat, 9UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 || mat(0,4) != 5 ||
          mat(0,5) != 6 || mat(0,6) != 7 || mat(0,7) != 8 || mat(0,8) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 3x3 initialization constructor
   {
      test_ = "Row-major StaticMatrix 3x3 initialization constructor";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 9UL );
      checkNonZeros( mat, 9UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ||
          mat(2,0) != 7 || mat(2,1) != 8 || mat(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n( 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 9x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 9x1 initialization constructor";

      blaze::StaticMatrix<int,9UL,1UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

      checkRows    ( mat, 9UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 9UL );
      checkNonZeros( mat, 9UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 || mat(4,0) != 5 ||
          mat(5,0) != 6 || mat(6,0) != 7 || mat(7,0) != 8 || mat(8,0) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n( 6 )\n( 7 )\n( 8 )\n( 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x10 initialization constructor
   {
      test_ = "Row-major StaticMatrix 1x10 initialization constructor";

      blaze::StaticMatrix<int,1UL,10UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );

      checkRows    ( mat,  1UL );
      checkColumns ( mat, 10UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat, 10UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 || mat(0,4) != 5 ||
          mat(0,5) != 6 || mat(0,6) != 7 || mat(0,7) != 8 || mat(0,8) != 9 || mat(0,9) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x5 initialization constructor
   {
      test_ = "Row-major StaticMatrix 2x5 initialization constructor";

      blaze::StaticMatrix<int,2UL,5UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat, 10UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 || mat(0,4) != 5 ||
          mat(1,0) != 6 || mat(1,1) != 7 || mat(1,2) != 8 || mat(1,3) != 9 || mat(1,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n( 6 7 8 9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 5x2 initialization constructor
   {
      test_ = "Row-major StaticMatrix 5x2 initialization constructor";

      blaze::StaticMatrix<int,5UL,2UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );

      checkRows    ( mat,  5UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat, 10UL );

      if( mat(0,0) != 1 || mat(0,1) !=  2 ||
          mat(1,0) != 3 || mat(1,1) !=  4 ||
          mat(2,0) != 5 || mat(2,1) !=  6 ||
          mat(3,0) != 7 || mat(3,1) !=  8 ||
          mat(4,0) != 9 || mat(4,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 )\n( 3 4 )\n( 5 6 )\n( 7 8 )\n( 9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 10x1 initialization constructor
   {
      test_ = "Row-major StaticMatrix 10x1 initialization constructor";

      blaze::StaticMatrix<int,10UL,1UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );

      checkRows    ( mat, 10UL );
      checkColumns ( mat,  1UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat, 10UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 || mat(4,0) != 5 ||
          mat(5,0) != 6 || mat(6,0) != 7 || mat(7,0) != 8 || mat(8,0) != 9 || mat(9,0) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  1 )\n(  2 )\n(  3 )\n(  4 )\n(  5 )\n"
                                     "(  6 )\n(  7 )\n(  8 )\n(  9 )\n( 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Copy constructor
   {
      test_ = "Row-major StaticMatrix copy constructor";

      blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat1( 1, 2, 3, 4, 5, 6 );
      blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );

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
   // Column-major matrix tests
   //=====================================================================================

   // Default constructor
   {
      test_ = "Column-major StaticMatrix default constructor";

      blaze::StaticMatrix<int,3UL,4UL,blaze::columnMajor> mat;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat,  0UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 || mat(1,3) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Homogeneous initialization
   {
      test_ = "Column-major StaticMatrix homogeneous initialization constructor";

      blaze::StaticMatrix<int,3UL,4UL,blaze::columnMajor> mat( 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );

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

   // 1x2 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x2 initialization constructor";

      blaze::StaticMatrix<int,1UL,2UL,blaze::columnMajor> mat( 1, 2 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 2x1 initialization constructor";

      blaze::StaticMatrix<int,2UL,1UL,blaze::columnMajor> mat( 1, 2 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x3 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x3 initialization constructor";

      blaze::StaticMatrix<int,1UL,3UL,blaze::columnMajor> mat( 1, 2, 3 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 3x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 3x1 initialization constructor";

      blaze::StaticMatrix<int,3UL,1UL,blaze::columnMajor> mat( 1, 2, 3 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x4 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x4 initialization constructor";

      blaze::StaticMatrix<int,1UL,4UL,blaze::columnMajor> mat( 1, 2, 3, 4 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x2 initialization constructor
   {
      test_ = "Column-major StaticMatrix 2x2 initialization constructor";

      blaze::StaticMatrix<int,2UL,2UL,blaze::columnMajor> mat( 1, 2, 3, 4 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( mat(0,0) != 1 || mat(0,1) != 3 || mat(1,0) != 2 || mat(1,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 3 )\n( 2 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 4x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 4x1 initialization constructor";

      blaze::StaticMatrix<int,4UL,1UL,blaze::columnMajor> mat( 1, 2, 3, 4 );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x5 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x5 initialization constructor";

      blaze::StaticMatrix<int,1UL,5UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 5UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 || mat(0,4) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 5x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 5x1 initialization constructor";

      blaze::StaticMatrix<int,5UL,1UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5 );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 5UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 || mat(4,0) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x6 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x6 initialization constructor";

      blaze::StaticMatrix<int,1UL,6UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 6UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(0,3) != 4 || mat(0,4) != 5 || mat(0,5) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x3 initialization constructor
   {
      test_ = "Column-major StaticMatrix 2x3 initialization constructor";

      blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

      if( mat(0,0) != 1 || mat(0,1) != 3 || mat(0,2) != 5 ||
          mat(1,0) != 2 || mat(1,1) != 4 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 3 5 )\n( 2 4 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

    // 3x2 initialization constructor
   {
      test_ = "Column-major StaticMatrix 3x2 initialization constructor";

      blaze::StaticMatrix<int,3UL,2UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

      if( mat(0,0) != 1 || mat(0,1) != 4 ||
          mat(1,0) != 2 || mat(1,1) != 5 ||
          mat(2,0) != 3 || mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 4 )\n( 2 5 )\n( 3 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 6x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 6x1 initialization constructor";

      blaze::StaticMatrix<int,6UL,1UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 ||
          mat(3,0) != 4 || mat(4,0) != 5 || mat(5,0) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n( 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x7 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x7 initialization constructor";

      blaze::StaticMatrix<int,1UL,7UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 7UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 ||
          mat(0,4) != 5 || mat(0,5) != 6 || mat(0,6) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 7x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 7x1 initialization constructor";

      blaze::StaticMatrix<int,7UL,1UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7 );

      checkRows    ( mat, 7UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 ||
          mat(4,0) != 5 || mat(5,0) != 6 || mat(6,0) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n( 6 )\n( 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x8 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x8 initialization constructor";

      blaze::StaticMatrix<int,1UL,8UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 8UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 8UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 ||
          mat(0,4) != 5 || mat(0,5) != 6 || mat(0,6) != 7 || mat(0,7) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x4 initialization constructor
   {
      test_ = "Column-major StaticMatrix 2x4 initialization constructor";

      blaze::StaticMatrix<int,2UL,4UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 8UL );

      if( mat(0,0) != 1 || mat(0,1) != 3 || mat(0,2) != 5 || mat(0,3) != 7 ||
          mat(1,0) != 2 || mat(1,1) != 4 || mat(1,2) != 6 || mat(1,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 3 5 7 )\n( 2 4 6 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 4x2 initialization constructor
   {
      test_ = "Column-major StaticMatrix 4x2 initialization constructor";

      blaze::StaticMatrix<int,4UL,2UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8 );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 8UL );

      if( mat(0,0) != 1 || mat(0,1) != 5 || mat(1,0) != 2 || mat(1,1) != 6 ||
          mat(2,0) != 3 || mat(2,1) != 7 || mat(3,0) != 4 || mat(3,1) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 5 )\n( 2 6 )\n( 3 7 )\n( 4 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 8x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 8x1 initialization constructor";

      blaze::StaticMatrix<int,8UL,1UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8 );

      checkRows    ( mat, 8UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 8UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 ||
          mat(4,0) != 5 || mat(5,0) != 6 || mat(6,0) != 7 || mat(7,0) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n( 6 )\n( 7 )\n( 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x9 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x9 initialization constructor";

      blaze::StaticMatrix<int,1UL,9UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

      checkRows    ( mat, 1UL );
      checkColumns ( mat, 9UL );
      checkCapacity( mat, 9UL );
      checkNonZeros( mat, 9UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 || mat(0,4) != 5 ||
          mat(0,5) != 6 || mat(0,6) != 7 || mat(0,7) != 8 || mat(0,8) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 3x3 initialization constructor
   {
      test_ = "Column-major StaticMatrix 3x3 initialization constructor";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 9UL );
      checkNonZeros( mat, 9UL );

      if( mat(0,0) != 1 || mat(0,1) != 4 || mat(0,2) != 7 ||
          mat(1,0) != 2 || mat(1,1) != 5 || mat(1,2) != 8 ||
          mat(2,0) != 3 || mat(2,1) != 6 || mat(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 4 7 )\n( 2 5 8 )\n( 3 6 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 9x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 9x1 initialization constructor";

      blaze::StaticMatrix<int,9UL,1UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

      checkRows    ( mat, 9UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 9UL );
      checkNonZeros( mat, 9UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 || mat(4,0) != 5 ||
          mat(5,0) != 6 || mat(6,0) != 7 || mat(7,0) != 8 || mat(8,0) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 )\n( 2 )\n( 3 )\n( 4 )\n( 5 )\n( 6 )\n( 7 )\n( 8 )\n( 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 1x10 initialization constructor
   {
      test_ = "Column-major StaticMatrix 1x10 initialization constructor";

      blaze::StaticMatrix<int,1UL,10UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );

      checkRows    ( mat,  1UL );
      checkColumns ( mat, 10UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat, 10UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 || mat(0,3) != 4 || mat(0,4) != 5 ||
          mat(0,5) != 6 || mat(0,6) != 7 || mat(0,7) != 8 || mat(0,8) != 9 || mat(0,9) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 2x5 initialization constructor
   {
      test_ = "Column-major StaticMatrix 2x5 initialization constructor";

      blaze::StaticMatrix<int,2UL,5UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat, 10UL );

      if( mat(0,0) != 1 || mat(0,1) != 3 || mat(0,2) != 5 || mat(0,3) != 7 || mat(0,4) != 9 ||
          mat(1,0) != 2 || mat(1,1) != 4 || mat(1,2) != 6 || mat(1,3) != 8 || mat(1,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 3 5 7 9 )\n( 2 4 6 8 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 5x2 initialization constructor
   {
      test_ = "Column-major StaticMatrix 5x2 initialization constructor";

      blaze::StaticMatrix<int,5UL,2UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );

      checkRows    ( mat,  5UL );
      checkColumns ( mat,  2UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat, 10UL );

      if( mat(0,0) != 1 || mat(0,1) !=  6 ||
          mat(1,0) != 2 || mat(1,1) !=  7 ||
          mat(2,0) != 3 || mat(2,1) !=  8 ||
          mat(3,0) != 4 || mat(3,1) !=  9 ||
          mat(4,0) != 5 || mat(4,1) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 6 )\n( 2 7 )\n( 3 8 )\n( 4 9 )\n( 5 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 10x1 initialization constructor
   {
      test_ = "Column-major StaticMatrix 10x1 initialization constructor";

      blaze::StaticMatrix<int,10UL,1UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );

      checkRows    ( mat, 10UL );
      checkColumns ( mat,  1UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat, 10UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,0) != 3 || mat(3,0) != 4 || mat(4,0) != 5 ||
          mat(5,0) != 6 || mat(6,0) != 7 || mat(7,0) != 8 || mat(8,0) != 9 || mat(9,0) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  1 )\n(  2 )\n(  3 )\n(  4 )\n(  5 )\n"
                                     "(  6 )\n(  7 )\n(  8 )\n(  9 )\n( 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Copy constructor
   {
      test_ = "Column-major StaticMatrix copy constructor";

      blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat1( 1, 2, 3, 4, 5, 6 );
      blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 3 || mat2(0,2) != 5 ||
          mat2(1,0) != 2 || mat2(1,1) != 4 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 3 5 )\n( 2 4 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the StaticMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void StaticMatrix::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major StaticMatrix::operator()";

      // Writing the first element
      blaze::StaticMatrix<int,3UL,5UL,blaze::rowMajor> mat;
      mat(2,1) = 1;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat,  1UL );

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
      test_ = "Column-major StaticMatrix::operator()";

      // Writing the first element
      blaze::StaticMatrix<int,3UL,5UL,blaze::columnMajor> mat;
      mat(2,1) = 1;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat,  1UL );

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
/*!\brief Test of the nonZeros member function of StaticMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of StaticMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticMatrix::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major StaticMatrix::nonZeros()";

      {
         blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat( 1, 2, 0, 3, 4, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 4UL );

         if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 0 ||
             mat(1,0) != 3 || mat(1,1) != 4 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 2 0 )\n( 3 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major StaticMatrix::nonZeros()";

      {
         blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat( 1, 2, 0, 3, 4, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 4UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 4 ||
             mat(1,0) != 2 || mat(1,1) != 3 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 4 )\n( 2 3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reset member function of StaticMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of StaticMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticMatrix::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major StaticMatrix::reset()";

      // Initialization check
      blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

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

      // Resetting the matrix
      mat.reset();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 0UL );

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
      test_ = "Column-major StaticMatrix::reset()";

      // Initialization check
      blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

      if( mat(0,0) != 1 || mat(0,1) != 3 || mat(0,2) != 5 ||
          mat(1,0) != 2 || mat(1,1) != 4 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the matrix
      mat.reset();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 0UL );

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
/*!\brief Test of the transpose member function of the StaticMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the transpose member functions of the StaticMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticMatrix::testTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major StaticMatrix::normalize()";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

      mat.transpose();

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 9UL );
      checkNonZeros( mat, 9UL );

      if( mat(0,0) != 1 || mat(0,1) != 4 || mat(0,2) != 7 ||
          mat(1,0) != 2 || mat(1,1) != 5 || mat(1,2) != 8 ||
          mat(2,0) != 3 || mat(2,1) != 6 || mat(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 4 7 )\n( 2 5 8 )\n( 3 6 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major StaticMatrix::normalize()";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

      mat.transpose();

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 9UL );
      checkNonZeros( mat, 9UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ||
          mat(2,0) != 7 || mat(2,1) != 8 || mat(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n( 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the isDiagonal member function of the StaticMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isDiagonal member functions of the StaticMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticMatrix::testIsDiagonal()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major StaticMatrix::isDiagonal()";

      // Non-quadratic matrix
      {
         blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( 1, 0, 0, 0, 2, 0, 0, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( 1, 0, 4, 0, 2, 0, 0, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );

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
      test_ = "Column-major StaticMatrix::isDiagonal()";

      // Non-quadratic matrix
      {
         blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( 1, 0, 0, 0, 2, 0, 0, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( 1, 0, 4, 0, 2, 0, 0, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );

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
/*!\brief Test of the isSymmetric member function of the StaticMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the isSymmetric member functions of the StaticMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticMatrix::testIsSymmetric()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major StaticMatrix::isSymmetric()";

      // Non-quadratic matrix
      {
         blaze::StaticMatrix<int,2UL,3UL,blaze::rowMajor> mat;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( 1, 0, 0, 0, 2, 0, 0, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( 1, 0, 4, 0, 2, 0, 0, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( 1, 0, 4, 0, 2, 0, 4, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );

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
      test_ = "Column-major StaticMatrix::isSymmetric()";

      // Non-quadratic matrix
      {
         blaze::StaticMatrix<int,2UL,3UL,blaze::columnMajor> mat;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( 1, 0, 0, 0, 2, 0, 0, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( 1, 0, 4, 0, 2, 0, 0, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );

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
         blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( 1, 0, 4, 0, 2, 0, 4, 0, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );

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
/*!\brief Test of the scale member function of StaticMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of StaticMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticMatrix::testScale()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major StaticMatrix::scale()";

      // Initialization check
      blaze::StaticMatrix<int,3UL,2UL,blaze::rowMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

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

      blaze::StaticMatrix<complex<float>,2UL,2UL,blaze::rowMajor> mat;
      mat(0,0) = complex<float>( 1.0F, 0.0F );
      mat(0,1) = complex<float>( 2.0F, 0.0F );
      mat(1,0) = complex<float>( 3.0F, 0.0F );
      mat(1,1) = complex<float>( 4.0F, 0.0F );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );

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
      test_ = "Column-major StaticMatrix::scale()";

      // Initialization check
      blaze::StaticMatrix<int,3UL,2UL,blaze::columnMajor> mat( 1, 2, 3, 4, 5, 6 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );

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

      blaze::StaticMatrix<complex<float>,2UL,2UL,blaze::columnMajor> mat;
      mat(0,0) = complex<float>( 1.0F, 0.0F );
      mat(0,1) = complex<float>( 2.0F, 0.0F );
      mat(1,0) = complex<float>( 3.0F, 0.0F );
      mat(1,1) = complex<float>( 4.0F, 0.0F );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );

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
/*!\brief Test of the swap functionality of the StaticMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the swap function of the StaticMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticMatrix::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major StaticMatrix swap";

      blaze::StaticMatrix<int,2UL,2UL,blaze::rowMajor> mat1( 1, 2, 0, 3 );
      blaze::StaticMatrix<int,2UL,2UL,blaze::rowMajor> mat2( 4, 3, 2, 1 );

      swap( mat1, mat2 );

      checkRows    ( mat1, 2UL );
      checkColumns ( mat1, 2UL );
      checkCapacity( mat1, 4UL );
      checkNonZeros( mat1, 4UL );

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
      test_ = "Column-major StaticMatrix swap";

      blaze::StaticMatrix<int,2UL,2UL,blaze::columnMajor> mat1( 1, 2, 0, 3 );
      blaze::StaticMatrix<int,2UL,2UL,blaze::columnMajor> mat2( 4, 3, 2, 1 );

      swap( mat1, mat2 );

      checkRows    ( mat1, 2UL );
      checkColumns ( mat1, 2UL );
      checkCapacity( mat1, 4UL );
      checkNonZeros( mat1, 4UL );

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

} // namespace staticmatrix

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
   std::cout << "   Running StaticMatrix test..." << std::endl;

   try
   {
      RUN_STATICMATRIX_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during StaticMatrix test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
