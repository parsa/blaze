//=================================================================================================
/*!
//  \file src/mathtest/lapack/OperationTest.cpp
//  \brief Source file for the LAPACK operation test
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
#include <sstream>
#include <stdexcept>
#include <blaze/math/LAPACK.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blazetest/mathtest/lapack/OperationTest.h>
#include <blazetest/system/LAPACK.h>


namespace blazetest {

namespace mathtest {

namespace lapack {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the OperationTest class test.
//
// \exception std::runtime_error Operation error detected.
*/
OperationTest::OperationTest()
{
   testLU();
   testQR();
   testInversion();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the LU factorization functionality.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the LU factorization functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void OperationTest::testLU()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   // Single precision matrices
   {
      test_ = "Row-major LU factorization (single precision)";

      blaze::StaticMatrix<float,3UL,3U,blaze::rowMajor> A( 2.0F, -1.0F, -2.0F,
                                                           4.0F,  1.0F, -7.0F,
                                                           6.0F,  3.0F, -8.0F );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::sgetrf( A, ipiv.data() );

      if( A(0,0) != 2.0F || A(0,1) != -0.5F || A(0,2) != -1.0F ||
          A(1,0) != 4.0F || A(1,1) !=  3.0F || A(1,2) != -1.0F ||
          A(2,0) != 6.0F || A(2,1) !=  6.0F || A(2,2) !=  4.0F ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( 2.0 -0.5 -1.0 )\n"
                                     "( 4.0  3.0 -1.0 )\n"
                                     "( 6.0  6.0  4.0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision matrices
   {
      test_ = "Row-major LU factorization (double precision)";

      blaze::StaticMatrix<double,3UL,3U,blaze::rowMajor> A( 2.0, -1.0, -2.0,
                                                            4.0,  1.0, -7.0,
                                                            6.0,  3.0, -8.0 );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::dgetrf( A, ipiv.data() );

      if( A(0,0) != 2.0 || A(0,1) != -0.5 || A(0,2) != -1.0 ||
          A(1,0) != 4.0 || A(1,1) !=  3.0 || A(1,2) != -1.0 ||
          A(2,0) != 6.0 || A(2,1) !=  6.0 || A(2,2) !=  4.0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( 2.0 -0.5 -1.0 )\n"
                                     "( 4.0  3.0 -1.0 )\n"
                                     "( 6.0  6.0  4.0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single precision complex matrices
   {
      test_ = "Row-major LU factorization (single precision complex)";

      typedef blaze::complex<float>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::rowMajor> A( cplx( 2.0F ), cplx( -1.0F ), cplx( -2.0F ),
                                                          cplx( 4.0F ), cplx(  1.0F ), cplx( -7.0F ),
                                                          cplx( 6.0F ), cplx(  3.0F ), cplx( -8.0F ) );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::cgetrf( A, ipiv.data() );

      if( A(0,0) != cplx( 2.0F ) || A(0,1) != cplx( -0.5F ) || A(0,2) != cplx( -1.0F ) ||
          A(1,0) != cplx( 4.0F ) || A(1,1) != cplx(  3.0F ) || A(1,2) != cplx( -1.0F ) ||
          A(2,0) != cplx( 6.0F ) || A(2,1) != cplx(  6.0F ) || A(2,2) != cplx(  4.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( (2.0,0.0) (-0.5,0.0) (-1.0,0.0) )\n"
                                     "( (4.0,0.0) ( 3.0,0.0) (-1.0,0.0) )\n"
                                     "( (6.0,0.0) ( 6.0,0.0) ( 4.0,0.0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision complex matrices
   {
      test_ = "Row-major LU factorization (double precision complex)";

      typedef blaze::complex<double>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::rowMajor> A( cplx( 2.0 ), cplx( -1.0 ), cplx( -2.0 ),
                                                          cplx( 4.0 ), cplx(  1.0 ), cplx( -7.0 ),
                                                          cplx( 6.0 ), cplx(  3.0 ), cplx( -8.0 ) );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::zgetrf( A, ipiv.data() );

      if( A(0,0) != cplx( 2.0 ) || A(0,1) != cplx( -0.5 ) || A(0,2) != cplx( -1.0 ) ||
          A(1,0) != cplx( 4.0 ) || A(1,1) != cplx(  3.0 ) || A(1,2) != cplx( -1.0 ) ||
          A(2,0) != cplx( 6.0 ) || A(2,1) != cplx(  6.0 ) || A(2,2) != cplx(  4.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( (2.0,0.0) (-0.5,0.0) (-1.0,0.0) )\n"
                                     "( (4.0,0.0) ( 3.0,0.0) (-1.0,0.0) )\n"
                                     "( (6.0,0.0) ( 6.0,0.0) ( 4.0,0.0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   // Single precision matrices
   {
      test_ = "Column-major LU factorization (single precision)";

      blaze::StaticMatrix<float,3UL,3U,blaze::columnMajor> A( 2.0F, -1.0F, -2.0F,
                                                              4.0F,  1.0F, -7.0F,
                                                              6.0F,  3.0F, -8.0F );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::sgetrf( A, ipiv.data() );

      if( A(0,0) !=  2.0F || A(0,1) !=  4.0F || A(0,2) != 6.0F ||
          A(1,0) != -0.5F || A(1,1) !=  3.0F || A(1,2) != 6.0F ||
          A(2,0) != -1.0F || A(2,1) != -1.0F || A(2,2) != 4.0F ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n(  2.0  4.0  6.0 )\n"
                                     "( -0.5  3.0  6.0 )\n"
                                     "( -1.0 -1.0  4.0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision matrices
   {
      test_ = "Column-major LU factorization (double precision)";

      blaze::StaticMatrix<double,3UL,3U,blaze::columnMajor> A( 2.0, -1.0, -2.0,
                                                               4.0,  1.0, -7.0,
                                                               6.0,  3.0, -8.0 );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::dgetrf( A, ipiv.data() );

      if( A(0,0) !=  2.0 || A(0,1) !=  4.0 || A(0,2) != 6.0 ||
          A(1,0) != -0.5 || A(1,1) !=  3.0 || A(1,2) != 6.0 ||
          A(2,0) != -1.0 || A(2,1) != -1.0 || A(2,2) != 4.0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n(  2.0  4.0  6.0 )\n"
                                     "( -0.5  3.0  6.0 )\n"
                                     "( -1.0 -1.0  4.0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single precision complex matrices
   {
      test_ = "Column-major LU factorization (single precision complex)";

      typedef blaze::complex<float>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::columnMajor> A( cplx( 2.0F ), cplx( -1.0F ), cplx( -2.0F ),
                                                             cplx( 4.0F ), cplx(  1.0F ), cplx( -7.0F ),
                                                             cplx( 6.0F ), cplx(  3.0F ), cplx( -8.0F ) );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::cgetrf( A, ipiv.data() );

      if( A(0,0) != cplx(  2.0F ) || A(0,1) != cplx(  4.0F ) || A(0,2) != cplx( 6.0F ) ||
          A(1,0) != cplx( -0.5F ) || A(1,1) != cplx(  3.0F ) || A(1,2) != cplx( 6.0F ) ||
          A(2,0) != cplx( -1.0F ) || A(2,1) != cplx( -1.0F ) || A(2,2) != cplx( 4.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( ( 2.0,0.0) ( 4.0,0.0) (6.0,0.0) )\n"
                                     "( (-0.5,0.0) ( 3.0,0.0) (6.0,0.0) )\n"
                                     "( (-1.0,0.0) (-1.0,0.0) (4.0,0.0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision complex matrices
   {
      test_ = "Column-major LU factorization (double precision complex)";

      typedef blaze::complex<double>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::columnMajor> A( cplx( 2.0 ), cplx( -1.0 ), cplx( -2.0 ),
                                                             cplx( 4.0 ), cplx(  1.0 ), cplx( -7.0 ),
                                                             cplx( 6.0 ), cplx(  3.0 ), cplx( -8.0 ) );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::zgetrf( A, ipiv.data() );

      if( A(0,0) != cplx(  2.0 ) || A(0,1) != cplx(  4.0 ) || A(0,2) != cplx( 6.0 ) ||
          A(1,0) != cplx( -0.5 ) || A(1,1) != cplx(  3.0 ) || A(1,2) != cplx( 6.0 ) ||
          A(2,0) != cplx( -1.0 ) || A(2,1) != cplx( -1.0 ) || A(2,2) != cplx( 4.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( ( 2.0,0.0) ( 4.0,0.0) (6.0,0.0) )\n"
                                     "( (-0.5,0.0) ( 3.0,0.0) (6.0,0.0) )\n"
                                     "( (-1.0,0.0) (-1.0,0.0) (4.0,0.0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the QR factorization functionality.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the QR factorization functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void OperationTest::testQR()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   // Single precision matrices
   {
      test_ = "Row-major QR factorization (single precision)";

      blaze::StaticMatrix<float,3UL,3U,blaze::rowMajor> A( 1.0F, 0.0F, 0.0F,
                                                           0.0F, 1.0F, 0.0F,
                                                           1.0F, 1.0F, 1.0F );
      blaze::StaticVector<float,3UL,blaze::columnVector> tau;

      blaze::sgeqrf( A, tau.data() );

      if( A(0,0) != 1.0F || A(0,1) != 0.0F || A(0,2) != 0.0F ||
          A(1,0) != 0.0F || A(1,1) != 1.0F || A(1,2) != 0.0F ||
          A(2,0) != 1.0F || A(2,1) != 1.0F || A(2,2) != 1.0F ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( 1  0  0 )\n"
                                     "( 0  1  0 )\n"
                                     "( 1  1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision matrices
   {
      test_ = "Row-major QR factorization (double precision)";

      blaze::StaticMatrix<double,3UL,3U,blaze::rowMajor> A( 1.0, 0.0, 0.0,
                                                            0.0, 1.0, 0.0,
                                                            1.0, 1.0, 1.0 );
      blaze::StaticVector<double,3UL,blaze::columnVector> tau;

      blaze::dgeqrf( A, tau.data() );

      if( A(0,0) != 1.0 || A(0,1) != 0.0 || A(0,2) != 0.0 ||
          A(1,0) != 0.0 || A(1,1) != 1.0 || A(1,2) != 0.0 ||
          A(2,0) != 1.0 || A(2,1) != 1.0 || A(2,2) != 1.0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( 1  0  0 )\n"
                                     "( 0  1  0 )\n"
                                     "( 1  1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single precision complex matrices
   {
      test_ = "Row-major QR factorization (single precision complex)";

      typedef blaze::complex<float>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::rowMajor> A( cplx( 1.0F ), cplx( 0.0F ), cplx( 0.0F ),
                                                          cplx( 0.0F ), cplx( 1.0F ), cplx( 0.0F ),
                                                          cplx( 1.0F ), cplx( 1.0F ), cplx( 1.0F ) );
      blaze::StaticVector<cplx,3UL,blaze::columnVector> tau;

      blaze::cgeqrf( A, tau.data() );

      if( A(0,0) != cplx( 1.0F ) || A(0,1) != cplx( 0.0F ) || A(0,2) != cplx( 0.0F ) ||
          A(1,0) != cplx( 0.0F ) || A(1,1) != cplx( 1.0F ) || A(1,2) != cplx( 0.0F ) ||
          A(2,0) != cplx( 1.0F ) || A(2,1) != cplx( 1.0F ) || A(2,2) != cplx( 1.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (1,0) (0,0) )\n"
                                     "( (1,0) (1,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision complex matrices
   {
      test_ = "Row-major QR factorization (double precision complex)";

      typedef blaze::complex<double>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::rowMajor> A( cplx( 1.0 ), cplx( 0.0 ), cplx( 0.0 ),
                                                          cplx( 0.0 ), cplx( 1.0 ), cplx( 0.0 ),
                                                          cplx( 1.0 ), cplx( 1.0 ), cplx( 1.0 ) );
      blaze::StaticVector<cplx,3UL,blaze::columnVector> tau;

      blaze::zgeqrf( A, tau.data() );

      if( A(0,0) != cplx( 1.0 ) || A(0,1) != cplx( 0.0 ) || A(0,2) != cplx( 0.0 ) ||
          A(1,0) != cplx( 0.0 ) || A(1,1) != cplx( 1.0 ) || A(1,2) != cplx( 0.0 ) ||
          A(2,0) != cplx( 1.0 ) || A(2,1) != cplx( 1.0 ) || A(2,2) != cplx( 1.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (1,0) (0,0) )\n"
                                     "( (1,0) (1,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   // Single precision matrices
   {
      test_ = "Row-major QR factorization (single precision)";

      blaze::StaticMatrix<float,3UL,3U,blaze::columnMajor> A( 1.0F, 0.0F, 0.0F,
                                                              0.0F, 1.0F, 0.0F,
                                                              1.0F, 1.0F, 1.0F );
      blaze::StaticVector<float,3UL,blaze::columnVector> tau;

      blaze::sgeqrf( A, tau.data() );

      if( A(0,0) != 1.0F || A(0,1) != 0.0F || A(0,2) != 1.0F ||
          A(1,0) != 0.0F || A(1,1) != 1.0F || A(1,2) != 1.0F ||
          A(2,0) != 0.0F || A(2,1) != 0.0F || A(2,2) != 1.0F ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( 1  0  1 )\n"
                                     "( 0  1  1 )\n"
                                     "( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision matrices
   {
      test_ = "Row-major QR factorization (double precision)";

      blaze::StaticMatrix<double,3UL,3U,blaze::columnMajor> A( 1.0, 0.0, 0.0,
                                                               0.0, 1.0, 0.0,
                                                               1.0, 1.0, 1.0 );
      blaze::StaticVector<double,3UL,blaze::columnVector> tau;

      blaze::dgeqrf( A, tau.data() );

      if( A(0,0) != 1.0 || A(0,1) != 0.0 || A(0,2) != 1.0 ||
          A(1,0) != 0.0 || A(1,1) != 1.0 || A(1,2) != 1.0 ||
          A(2,0) != 0.0 || A(2,1) != 0.0 || A(2,2) != 1.0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( 1  0  1 )\n"
                                     "( 0  1  1 )\n"
                                     "( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single precision complex matrices
   {
      test_ = "Row-major QR factorization (single precision complex)";

      typedef blaze::complex<float>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::columnMajor> A( cplx( 1.0F ), cplx( 0.0F ), cplx( 0.0F ),
                                                             cplx( 0.0F ), cplx( 1.0F ), cplx( 0.0F ),
                                                             cplx( 1.0F ), cplx( 1.0F ), cplx( 1.0F ) );
      blaze::StaticVector<cplx,3UL,blaze::columnVector> tau;

      blaze::cgeqrf( A, tau.data() );

      if( A(0,0) != cplx( 1.0F ) || A(0,1) != cplx( 0.0F ) || A(0,2) != cplx( 1.0F ) ||
          A(1,0) != cplx( 0.0F ) || A(1,1) != cplx( 1.0F ) || A(1,2) != cplx( 1.0F ) ||
          A(2,0) != cplx( 0.0F ) || A(2,1) != cplx( 0.0F ) || A(2,2) != cplx( 1.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( (1,0) (0,0) (1,0) )\n"
                                     "( (0,0) (1,0) (1,0) )\n"
                                     "( (0,0) (0,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision complex matrices
   {
      test_ = "Row-major QR factorization (double precision complex)";

      typedef blaze::complex<double>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::columnMajor> A( cplx( 1.0 ), cplx( 0.0 ), cplx( 0.0 ),
                                                             cplx( 0.0 ), cplx( 1.0 ), cplx( 0.0 ),
                                                             cplx( 1.0 ), cplx( 1.0 ), cplx( 1.0 ) );
      blaze::StaticVector<cplx,3UL,blaze::columnVector> tau;

      blaze::zgeqrf( A, tau.data() );

      if( A(0,0) != cplx( 1.0 ) || A(0,1) != cplx( 0.0 ) || A(0,2) != cplx( 1.0 ) ||
          A(1,0) != cplx( 0.0 ) || A(1,1) != cplx( 1.0 ) || A(1,2) != cplx( 1.0 ) ||
          A(2,0) != cplx( 0.0 ) || A(2,1) != cplx( 0.0 ) || A(2,2) != cplx( 1.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR factorization failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( (1,0) (0,0) (1,0) )\n"
                                     "( (0,0) (1,0) (1,0) )\n"
                                     "( (0,0) (0,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the matrix inversion functionality.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the matrix inversion functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void OperationTest::testInversion()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   // Single precision matrices
   {
      test_ = "Row-major inversion (single precision)";

      blaze::StaticMatrix<float,3UL,3U,blaze::rowMajor> A( 1.0F, 0.0F, 0.0F,
                                                           0.0F, 1.0F, 0.0F,
                                                           1.0F, 1.0F, 1.0F );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::sgetrf( A, ipiv.data() );
      blaze::sgetri( A, ipiv.data() );

      if( A(0,0) !=  1.0F || A(0,1) !=  0.0F || A(0,2) != 0.0F ||
          A(1,0) !=  0.0F || A(1,1) !=  1.0F || A(1,2) != 0.0F ||
          A(2,0) != -1.0F || A(2,1) != -1.0F || A(2,2) != 1.0F ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n(  1  0  0 )\n"
                                     "(  0  1  0 )\n"
                                     "( -1 -1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision matrices
   {
      test_ = "Row-major inversion (double precision)";

      blaze::StaticMatrix<double,3UL,3U,blaze::rowMajor> A( 1.0, 0.0, 0.0,
                                                            0.0, 1.0, 0.0,
                                                            1.0, 1.0, 1.0 );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::dgetrf( A, ipiv.data() );
      blaze::dgetri( A, ipiv.data() );

      if( A(0,0) !=  1.0 || A(0,1) !=  0.0 || A(0,2) != 0.0 ||
          A(1,0) !=  0.0 || A(1,1) !=  1.0 || A(1,2) != 0.0 ||
          A(2,0) != -1.0 || A(2,1) != -1.0 || A(2,2) != 1.0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n(  1  0  0 )\n"
                                     "(  0  1  0 )\n"
                                     "( -1 -1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single precision complex matrices
   {
      test_ = "Row-major matrix inversion (single precision complex)";

      typedef blaze::complex<float>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::rowMajor> A( cplx( 1.0F ), cplx( 0.0F ), cplx( 0.0F ),
                                                          cplx( 0.0F ), cplx( 1.0F ), cplx( 0.0F ),
                                                          cplx( 1.0F ), cplx( 1.0F ), cplx( 1.0F ) );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::cgetrf( A, ipiv.data() );
      blaze::cgetri( A, ipiv.data() );

      if( A(0,0) != cplx(  1.0F ) || A(0,1) != cplx(  0.0F ) || A(0,2) != cplx( 0.0F ) ||
          A(1,0) != cplx(  0.0F ) || A(1,1) != cplx(  1.0F ) || A(1,2) != cplx( 0.0F ) ||
          A(2,0) != cplx( -1.0F ) || A(2,1) != cplx( -1.0F ) || A(2,2) != cplx( 1.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( ( 1,0) ( 0,0) (0,0) )\n"
                                     "( ( 0,0) ( 3,0) (0,0) )\n"
                                     "( (-1,0) (-1,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision complex matrices
   {
      test_ = "Row-major matrix inversion (double precision complex)";

      typedef blaze::complex<double>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::rowMajor> A( cplx( 1.0 ), cplx( 0.0 ), cplx( 0.0 ),
                                                          cplx( 0.0 ), cplx( 1.0 ), cplx( 0.0 ),
                                                          cplx( 1.0 ), cplx( 1.0 ), cplx( 1.0 ) );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::zgetrf( A, ipiv.data() );
      blaze::zgetri( A, ipiv.data() );

      if( A(0,0) != cplx(  1.0 ) || A(0,1) != cplx(  0.0 ) || A(0,2) != cplx( 0.0 ) ||
          A(1,0) != cplx(  0.0 ) || A(1,1) != cplx(  1.0 ) || A(1,2) != cplx( 0.0 ) ||
          A(2,0) != cplx( -1.0 ) || A(2,1) != cplx( -1.0 ) || A(2,2) != cplx( 1.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( ( 1,0) ( 0,0) (0,0) )\n"
                                     "( ( 0,0) ( 3,0) (0,0) )\n"
                                     "( (-1,0) (-1,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   // Single precision matrices
   {
      test_ = "Column-major inversion (single precision)";

      blaze::StaticMatrix<float,3UL,3U,blaze::columnMajor> A( 1.0F, 0.0F, 1.0F,
                                                              0.0F, 1.0F, 1.0F,
                                                              0.0F, 0.0F, 1.0F );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::sgetrf( A, ipiv.data() );
      blaze::sgetri( A, ipiv.data() );

      if( A(0,0) !=  1.0F || A(0,1) !=  0.0F || A(0,2) != 0.0F ||
          A(1,0) !=  0.0F || A(1,1) !=  1.0F || A(1,2) != 0.0F ||
          A(2,0) != -1.0F || A(2,1) != -1.0F || A(2,2) != 1.0F ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n(  1  0  0 )\n"
                                     "(  0  1  0 )\n"
                                     "( -1 -1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision matrices
   {
      test_ = "Column-major inversion (double precision)";

      blaze::StaticMatrix<double,3UL,3U,blaze::columnMajor> A( 1.0, 0.0, 1.0,
                                                               0.0, 1.0, 1.0,
                                                               0.0, 0.0, 1.0 );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::dgetrf( A, ipiv.data() );
      blaze::dgetri( A, ipiv.data() );

      if( A(0,0) !=  1.0 || A(0,1) !=  0.0 || A(0,2) != 0.0 ||
          A(1,0) !=  0.0 || A(1,1) !=  1.0 || A(1,2) != 0.0 ||
          A(2,0) != -1.0 || A(2,1) != -1.0 || A(2,2) != 1.0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n(  1  0  0 )\n"
                                     "(  0  1  0 )\n"
                                     "( -1 -1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single precision complex matrices
   {
      test_ = "Column-major matrix inversion (single precision complex)";

      typedef blaze::complex<float>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::columnMajor> A( cplx( 1.0F ), cplx( 0.0F ), cplx( 1.0F ),
                                                             cplx( 0.0F ), cplx( 1.0F ), cplx( 1.0F ),
                                                             cplx( 0.0F ), cplx( 0.0F ), cplx( 1.0F ) );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::cgetrf( A, ipiv.data() );
      blaze::cgetri( A, ipiv.data() );

      if( A(0,0) != cplx(  1.0F ) || A(0,1) != cplx(  0.0F ) || A(0,2) != cplx( 0.0F ) ||
          A(1,0) != cplx(  0.0F ) || A(1,1) != cplx(  1.0F ) || A(1,2) != cplx( 0.0F ) ||
          A(2,0) != cplx( -1.0F ) || A(2,1) != cplx( -1.0F ) || A(2,2) != cplx( 1.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( ( 1,0) ( 0,0) (0,0) )\n"
                                     "( ( 0,0) ( 3,0) (0,0) )\n"
                                     "( (-1,0) (-1,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Double precision complex matrices
   {
      test_ = "Column-major matrix inversion (double precision complex)";

      typedef blaze::complex<double>  cplx;

      blaze::StaticMatrix<cplx,3UL,3U,blaze::columnMajor> A( cplx( 1.0 ), cplx( 0.0 ), cplx( 1.0 ),
                                                             cplx( 0.0 ), cplx( 1.0 ), cplx( 1.0 ),
                                                             cplx( 0.0 ), cplx( 0.0 ), cplx( 1.0 ) );
      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      blaze::zgetrf( A, ipiv.data() );
      blaze::zgetri( A, ipiv.data() );

      if( A(0,0) != cplx(  1.0 ) || A(0,1) != cplx(  0.0 ) || A(0,2) != cplx( 0.0 ) ||
          A(1,0) != cplx(  0.0 ) || A(1,1) != cplx(  1.0 ) || A(1,2) != cplx( 0.0 ) ||
          A(2,0) != cplx( -1.0 ) || A(2,1) != cplx( -1.0 ) || A(2,2) != cplx( 1.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Result:\n" << A << "\n"
             << "   Expected result:\n( ( 1,0) ( 0,0) (0,0) )\n"
                                     "( ( 0,0) ( 3,0) (0,0) )\n"
                                     "( (-1,0) (-1,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************

} // namespace lapack

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
   std::cout << "   Running LAPACK operation test..." << std::endl;

   try
   {
      RUN_LAPACK_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during LAPACK operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
