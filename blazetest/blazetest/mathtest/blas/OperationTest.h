//=================================================================================================
/*!
//  \file blazetest/mathtest/blas/OperationTest.h
//  \brief Header file for the BLAS operation test
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

#ifndef _BLAZETEST_MATHTEST_BLAS_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_BLAS_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/BLAS.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/UpperMatrix.h>
#include <blazetest/system/BLAS.h>


namespace blazetest {

namespace mathtest {

namespace blas {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the BLAS functionality.
//
// This class represents a test suite for BLAS functionality wrapped by the Blaze library.
*/
class OperationTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   template< typename Type > void testTrsm();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the LU decomposition functions (getrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the LU decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testTrsm()
{
#if BLAZETEST_MATHTEST_BLAS_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major triangular LSE (single right-hand side, left side, lower part)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsm( A, x, CblasLeft, CblasLower, 1.0 );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, left side, upper part)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsm( A, x, CblasLeft, CblasUpper, 1.0 );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, right side, lower part)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> b, x;
      randomize( b );

      x = b;

      blaze::trsm( A, x, CblasRight, CblasLower, 1.0 );

      if( ( x * A ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   x * A:\n" << ( x * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, right side, upper part)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> b, x;
      randomize( b );

      x = b;

      blaze::trsm( A, x, CblasRight, CblasUpper, 1.0 );

      if( ( x * A ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   x * A:\n" << ( x * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (multiple right-hand sides, left side, lower part)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnVector> B, X;
      randomize( B );

      X = B;

      blaze::trsm( A, X, CblasLeft, CblasLower, 1.0 );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (multiple right-hand sides, left side, upper part)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnVector> B, X;
      randomize( B );

      X = B;

      blaze::trsm( A, X, CblasLeft, CblasUpper, 1.0 );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (multiple right-hand sides, right side, lower part)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::columnVector> B, X;
      randomize( B );

      X = B;

      blaze::trsm( A, X, CblasRight, CblasLower, 1.0 );

      if( ( X * A ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   X * A:\n" << ( X * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (multiple right-hand sides, right side, upper part)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::columnVector> B, X;
      randomize( B );

      X = B;

      blaze::trsm( A, X, CblasRight, CblasUpper, 1.0 );

      if( ( X * A ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   X * A:\n" << ( X * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major triangular LSE (single right-hand side, left side, lower part)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsm( A, x, CblasLeft, CblasLower, 1.0 );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (single right-hand side, left side, upper part)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsm( A, x, CblasLeft, CblasUpper, 1.0 );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (single right-hand side, right side, lower part)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> b, x;
      randomize( b );

      x = b;

      blaze::trsm( A, x, CblasRight, CblasLower, 1.0 );

      if( ( x * A ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   x * A:\n" << ( x * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (single right-hand side, right side, upper part)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> b, x;
      randomize( b );

      x = b;

      blaze::trsm( A, x, CblasRight, CblasUpper, 1.0 );

      if( ( x * A ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   x * A:\n" << ( x * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (multiple right-hand sides, left side, lower part)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnVector> B, X;
      randomize( B );

      X = B;

      blaze::trsm( A, X, CblasLeft, CblasLower, 1.0 );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (multiple right-hand sides, left side, upper part)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnVector> B, X;
      randomize( B );

      X = B;

      blaze::trsm( A, X, CblasLeft, CblasUpper, 1.0 );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (multiple right-hand sides, right side, lower part)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::columnVector> B, X;
      randomize( B );

      X = B;

      blaze::trsm( A, X, CblasRight, CblasLower, 1.0 );

      if( ( X * A ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   X * A:\n" << ( X * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (multiple right-hand sides, right side, upper part)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::columnVector> B, X;
      randomize( B );

      X = B;

      blaze::trsm( A, X, CblasRight, CblasUpper, 1.0 );

      if( ( X * A ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   X * A:\n" << ( X * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the BLAS functionality.
//
// \return void
*/
void runTest()
{
   OperationTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the BLAS operation test.
*/
#define RUN_BLAS_OPERATION_TEST \
   blazetest::mathtest::blas::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace blas

} // namespace mathtest

} // namespace blazetest

#endif
