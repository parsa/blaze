//=================================================================================================
/*!
//  \file blazetest/mathtest/lapack/OperationTest.h
//  \brief Header file for the LAPACK operation test
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

#ifndef _BLAZETEST_MATHTEST_LAPACK_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_LAPACK_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LAPACK.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blazetest/system/LAPACK.h>


namespace blazetest {

namespace mathtest {

namespace lapack {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the LAPACK functionality.
//
// This class represents a test suite for LAPACK functionality wrapped by the Blaze library.
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
   template< typename Type > void testGetrf();
   template< typename Type > void testSytrf();
   template< typename Type > void testHetrf();
   template< typename Type > void testPotrf();

   template< typename Type > void testGetri();
   template< typename Type > void testSytri();
   template< typename Type > void testHetri();
   template< typename Type > void testPotri();
   template< typename Type > void testTrtri();

   template< typename Type > void testGetrs();
   template< typename Type > void testSytrs();
   template< typename Type > void testHetrs();
   template< typename Type > void testPotrs();
   template< typename Type > void testTrtrs();

   template< typename Type > void testGesv();
   template< typename Type > void testSysv();
   template< typename Type > void testHesv();
   template< typename Type > void testPosv();
   template< typename Type > void testTrsv();

   template< typename Type > void testGeqrf();
   template< typename Type > void testOrgqr();
   template< typename Type > void testUngqr();

   template< typename Type > void testGerqf();
   template< typename Type > void testOrgrq();
   template< typename Type > void testUngrq();

   template< typename Type > void testGeqlf();
   template< typename Type > void testOrgql();
   template< typename Type > void testUngql();

   template< typename Type > void testGelqf();
   template< typename Type > void testOrglq();
   template< typename Type > void testUnglq();
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
void OperationTest::testGetrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "LU decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<int,2UL,blaze::columnVector> ipivA;
      blaze::StaticVector<int,2UL,blaze::columnVector> ipivB;

      blaze::getrf( A, ipivA.data() );
      blaze::getrf( B, ipivB.data() );

      if( A != trans( B ) || ipivA != ipivB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major pivot elements:\n" << ipivA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major pivot elements:\n" << ipivB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<int,2UL,blaze::columnVector> ipivA;
      blaze::StaticVector<int,2UL,blaze::columnVector> ipivB;

      blaze::getrf( A, ipivA.data() );
      blaze::getrf( B, ipivB.data() );

      if( A != trans( B ) || ipivA != ipivB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major pivot elements:\n" << ipivA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major pivot elements:\n" << ipivB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Bunch-Kaufman decomposition functions for symmetric matrices (sytrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Bunch-Kaufman decomposition functions for symmetric
// indefinite matrices for various data types. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >
void OperationTest::testSytrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Symmetric matrix decomposition";

   {
      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipivA;
      blaze::StaticVector<int,3UL,blaze::rowVector> ipivB;

      blaze::sytrf( A, 'L', ipivA.data() );
      blaze::sytrf( B, 'U', ipivB.data() );

      if( A != trans( B ) || ipivA != ipivB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major pivot elements:\n" << ipivA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major pivot elements:\n" << ipivB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Bunch-Kaufman decomposition functions for Hermitian matrices (hetrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Bunch-Kaufman decomposition functions for Hermitian
// indefinite matrices for various data types. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >
void OperationTest::testHetrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Hermitian matrix decomposition";

   {
      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipivA;
      blaze::StaticVector<int,3UL,blaze::rowVector> ipivB;

      blaze::hetrf( A, 'L', ipivA.data() );
      blaze::hetrf( B, 'U', ipivB.data() );

      if( A != ctrans( B ) || ipivA != ipivB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major pivot elements:\n" << ipivA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major pivot elements:\n" << ipivB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Cholesky decomposition functions (potrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Cholesky decomposition functions for various data types.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testPotrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Cholesky decomposition";

   {
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::potrf( A, 'L' );
      blaze::potrf( B, 'L' );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Column-major decomposition:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::potrf( A, 'U' );
      blaze::potrf( B, 'U' );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Column-major decomposition:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LU-based matrix inversion functions (getri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the LU-based matrix inversion functions for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testGetri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LU-based matrix inversion";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> Ainv( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::getrf( Ainv, ipiv.data() );
      blaze::getri( Ainv, ipiv.data() );

      if( !blaze::isIdentity( Ainv * A ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << Ainv << "\n"
             << "   Ainv * A = " << ( Ainv * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LU-based matrix inversion";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> Ainv( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::getrf( Ainv, ipiv.data() );
      blaze::getri( Ainv, ipiv.data() );

      if( !blaze::isIdentity( Ainv * A ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << Ainv << "\n"
             << "   Ainv * A = " << ( Ainv * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Bunch-Kaufman-based matrix inversion functions for symmetric matrices (sytri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Bunch-Kaufman-based matrix inversion functions for
// symmetric indefinite matrices for various data types. In case an error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testSytri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major symmetric matrix inversion (lower part)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::sytrf( B, 'L', ipiv.data() );
      blaze::sytri( B, 'L', ipiv.data() );

      B(0,1) = B(1,0);
      B(0,2) = B(2,0);
      B(1,2) = B(2,1);

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric matrix inversion (upper part)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::sytrf( B, 'U', ipiv.data() );
      blaze::sytri( B, 'U', ipiv.data() );

      B(1,0) = B(0,1);
      B(2,0) = B(0,2);
      B(2,1) = B(1,2);

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major symmetric matrix inversion (lower part)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::sytrf( B, 'L', ipiv.data() );
      blaze::sytri( B, 'L', ipiv.data() );

      B(0,1) = B(1,0);
      B(0,2) = B(2,0);
      B(1,2) = B(2,1);

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major symmetric matrix inversion (upper part)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::sytrf( B, 'U', ipiv.data() );
      blaze::sytri( B, 'U', ipiv.data() );

      B(1,0) = B(0,1);
      B(2,0) = B(0,2);
      B(2,1) = B(1,2);

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Bunch-Kaufman-based matrix inversion functions for Hermitian matrices (hetri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Bunch-Kaufman-based matrix inversion functions for
// Hermitian indefinite matrices for various data types. In case an error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testHetri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::conj;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Hermitian matrix inversion (lower part)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::hetrf( B, 'L', ipiv.data() );
      blaze::hetri( B, 'L', ipiv.data() );

      B(0,1) = conj( B(1,0) );
      B(0,2) = conj( B(2,0) );
      B(1,2) = conj( B(2,1) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian matrix inversion (upper part)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::hetrf( B, 'U', ipiv.data() );
      blaze::hetri( B, 'U', ipiv.data() );

      B(1,0) = conj( B(0,1) );
      B(2,0) = conj( B(0,2) );
      B(2,1) = conj( B(1,2) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Hermitian matrix inversion (lower part)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::hetrf( B, 'L', ipiv.data() );
      blaze::hetri( B, 'L', ipiv.data() );

      B(0,1) = conj( B(1,0) );
      B(0,2) = conj( B(2,0) );
      B(1,2) = conj( B(2,1) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Hermitian matrix inversion (upper part)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );
      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      blaze::hetrf( B, 'U', ipiv.data() );
      blaze::hetri( B, 'U', ipiv.data() );

      B(1,0) = conj( B(0,1) );
      B(2,0) = conj( B(0,2) );
      B(2,1) = conj( B(1,2) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Cholesky-based matrix inversion functions (potri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Cholesky-based matrix inversion functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testPotri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::conj;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Cholesky-based matrix inversion (lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::potrf( A, 'L' );
      blaze::potri( A, 'L' );

      A(0,1) = conj( A(1,0) );
      A(0,2) = conj( A(2,0) );
      A(1,2) = conj( A(2,1) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Cholesky-based matrix inversion (upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::potrf( A, 'U' );
      blaze::potri( A, 'U' );

      A(1,0) = conj( A(0,1) );
      A(2,0) = conj( A(0,2) );
      A(2,1) = conj( A(1,2) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Cholesky-based matrix inversion (lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::potrf( A, 'L' );
      blaze::potri( A, 'L' );

      A(0,1) = conj( A(1,0) );
      A(0,2) = conj( A(2,0) );
      A(1,2) = conj( A(2,1) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Cholesky-based matrix inversion (upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::potrf( A, 'U' );
      blaze::potri( A, 'U' );

      A(1,0) = conj( A(0,1) );
      A(2,0) = conj( A(0,2) );
      A(2,1) = conj( A(1,2) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the triangular matrix inversion functions (trtri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the triangular matrix inversion functions for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testTrtri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major lower triangular matrix inversion";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::trtri( B, 'L', 'N' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower triangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major lower unitriangular matrix inversion";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::trtri( B, 'L', 'U' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower unitriangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major upper triangular matrix inversion";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::trtri( B, 'U', 'N' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper triangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major upper unitriangular matrix inversion";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::trtri( B, 'U', 'U' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper unitriangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major lower triangular matrix inversion";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::trtri( B, 'L', 'N' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower triangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major lower unitriangular matrix inversion";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::trtri( B, 'L', 'U' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower unitriangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major upper triangular matrix inversion";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::trtri( B, 'U', 'N' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper triangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major upper unitriangular matrix inversion";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::trtri( B, 'U', 'U' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper unitriangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the general substitution functions (getrs).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the general substitution functions for various data types.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testGetrs()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major general LSE substitution (single right-hand side, not transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, x, 'N', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general LSE substitution (single right-hand side, transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = trans( A );
      x = b;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, x, 'T', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general LSE substitution (single right-hand side, conjugate transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = ctrans( A );
      x = b;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, x, 'C', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general LSE substitution (multiple right-hand sides, not transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, X, 'N', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general LSE substitution (multiple right-hand sides, transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = trans( A );
      X = B;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, X, 'T', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general LSE substitution (multiple right-hand sides, conjugate transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = ctrans( A );
      X = B;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, X, 'C', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major general LSE substitution (single right-hand side, not transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, x, 'N', ipiv.data() );

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
      test_ = "Column-major general LSE substitution (single right-hand side, transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = trans( A );
      x = b;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, x, 'T', ipiv.data() );

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
      test_ = "Column-major general LSE substitution (single right-hand side, conjugate transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = ctrans( A );
      x = b;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, x, 'C', ipiv.data() );

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
      test_ = "Column-major general LSE substitution (multiple right-hand sides, not transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, X, 'N', ipiv.data() );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * x:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major general LSE substitution (multiple right-hand sides, transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = trans( A );
      X = B;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, X, 'T', ipiv.data() );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * x:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major general LSE substitution (multiple right-hand sides, conjugate transposed)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = ctrans( A );
      X = B;

      blaze::getrf( LU, ipiv.data() );
      blaze::getrs( LU, X, 'C', ipiv.data() );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * x:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the symmetric indefinite substitution functions (sytrs).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the symmetric indefinite substitution functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testSytrs()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major symmetric indefinite LSE substitution (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::sytrf( LU, 'L', ipiv.data() );
      blaze::sytrs( LU, x, 'L', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE substitution (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::sytrf( LU, 'U', ipiv.data() );
      blaze::sytrs( LU, x, 'U', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE substitution (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::sytrf( LU, 'L', ipiv.data() );
      blaze::sytrs( LU, X, 'L', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE substitution (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::sytrf( LU, 'U', ipiv.data() );
      blaze::sytrs( LU, X, 'U', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major symmetric indefinite LSE substitution (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::sytrf( LU, 'L', ipiv.data() );
      blaze::sytrs( LU, x, 'L', ipiv.data() );

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
      test_ = "Column-major symmetric indefinite LSE substitution (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::sytrf( LU, 'U', ipiv.data() );
      blaze::sytrs( LU, x, 'U', ipiv.data() );

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
      test_ = "Column-major symmetric indefinite LSE substitution (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::sytrf( LU, 'L', ipiv.data() );
      blaze::sytrs( LU, X, 'L', ipiv.data() );

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
      test_ = "Column-major symmetric indefinite LSE substitution (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::sytrf( LU, 'U', ipiv.data() );
      blaze::sytrs( LU, X, 'U', ipiv.data() );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Hermitian indefinite substitution functions (hetrs).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Hermitian indefinite substitution functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testHetrs()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Hermitian indefinite LSE substitution (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::hetrf( LU, 'L', ipiv.data() );
      blaze::hetrs( LU, x, 'L', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE substitution (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::hetrf( LU, 'U', ipiv.data() );
      blaze::hetrs( LU, x, 'U', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE substitution (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::hetrf( LU, 'L', ipiv.data() );
      blaze::hetrs( LU, X, 'L', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE substitution (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::hetrf( LU, 'U', ipiv.data() );
      blaze::hetrs( LU, X, 'U', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Hermitian indefinite LSE substitution (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::hetrf( LU, 'L', ipiv.data() );
      blaze::hetrs( LU, x, 'L', ipiv.data() );

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
      test_ = "Column-major Hermitian indefinite LSE substitution (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      x = b;

      blaze::hetrf( LU, 'U', ipiv.data() );
      blaze::hetrs( LU, x, 'U', ipiv.data() );

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
      test_ = "Column-major Hermitian indefinite LSE substitution (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::hetrf( LU, 'L', ipiv.data() );
      blaze::hetrs( LU, X, 'L', ipiv.data() );

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
      test_ = "Column-major Hermitian indefinite LSE substitution (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::rowVector> ipiv;

      LU = A;
      X = B;

      blaze::hetrf( LU, 'U', ipiv.data() );
      blaze::hetrs( LU, X, 'U', ipiv.data() );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the positive definite substitution functions (potrs).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the positive definite substitution functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testPotrs()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major positive definite LSE substitution (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::potrf( LU, 'L' );
      blaze::potrs( LU, x, 'L' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE substitution (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::potrf( LU, 'U' );
      blaze::potrs( LU, x, 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE substitution (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::potrf( LU, 'L' );
      blaze::potrs( LU, X, 'L' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE substitution (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::potrf( LU, 'U' );
      blaze::potrs( LU, X, 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major positive definite LSE substitution (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::potrf( LU, 'L' );
      blaze::potrs( LU, x, 'L' );

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
      test_ = "Column-major positive definite LSE substitution (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::potrf( LU, 'U' );
      blaze::potrs( LU, x, 'U' );

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
      test_ = "Column-major positive definite LSE substitution (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::potrf( LU, 'L' );
      blaze::potrs( LU, X, 'L' );

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
      test_ = "Column-major positive definite LSE substitution (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::potrf( LU, 'U' );
      blaze::potrs( LU, X, 'U' );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the triangular substitution functions (trtrs).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the triangular substitution functions for various data types.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testTrtrs()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major triangular LSE substitution (single right-hand side, lower part, not transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'N', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (single right-hand side, lower part, transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'T', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (single right-hand side, lower part, conjugate transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'C', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (single right-hand side, lower part, not transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'N', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (single right-hand side, lower part, transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'T', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (single right-hand side, lower part, conjugate transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'C', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (single right-hand side, upper part, not transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'N', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (single right-hand side, upper part, transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'T', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (single right-hand side, upper part, conjugate transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'C', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (single right-hand side, upper part, not transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'N', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (single right-hand side, upper part, transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'T', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (single right-hand side, upper part, conjugate transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'C', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (multiple right-hand sides, lower part, not transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'N', 'N' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (multiple right-hand sides, lower part, transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'T', 'N' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (multiple right-hand sides, lower part, conjugate transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'C', 'N' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (multiple right-hand sides, lower part, not transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'N', 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (multiple right-hand sides, lower part, transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'T', 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (multiple right-hand sides, lower part, conjugate transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'C', 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (multiple right-hand sides, upper part, not transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'N', 'N' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (multiple right-hand sides, upper part, transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      LU = trans( A );
      X = B;

      blaze::trtrs( LU, X, 'L', 'T', 'N' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE substitution (multiple right-hand sides, upper part, conjugate transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      LU = ctrans( A );
      X = B;

      blaze::trtrs( LU, X, 'L', 'C', 'N' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (multiple right-hand sides, upper part, not transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'N', 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (multiple right-hand sides, upper part, transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'T', 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE substitution (multiple right-hand sides, upper part, conjugate transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticMatrix<Type,6UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'C', 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major triangular LSE substitution (single right-hand side, lower part, not transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'N', 'N' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE substitution (single right-hand side, lower part, transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'T', 'N' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE substitution (single right-hand side, lower part, conjugate transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'C', 'N' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE substitution (single right-hand side, lower part, not transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'N', 'U' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE substitution (single right-hand side, lower part, transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'T', 'U' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE substitution (single right-hand side, lower part, conjugate transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'C', 'U' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE substitution (single right-hand side, upper part, not transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'N', 'N' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE substitution (single right-hand side, upper part, transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'T', 'N' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE substitution (single right-hand side, upper part, conjugate transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'C', 'N' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE substitution (single right-hand side, upper part, not transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'U', 'N', 'U' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE substitution (single right-hand side, upper part, transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'T', 'U' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE substitution (single right-hand side, upper part, conjugate transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trtrs( LU, x, 'L', 'C', 'U' );

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
             << "   trans( A ) * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE substitution (multiple right-hand sides, lower part, not transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'N', 'N' );

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
      test_ = "Column-major triangular LSE substitution (multiple right-hand sides, lower part, transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'T', 'N' );

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
      test_ = "Column-major triangular LSE substitution (multiple right-hand sides, lower part, conjugate transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'C', 'N' );

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
      test_ = "Column-major unitriangular LSE substitution (multiple right-hand sides, lower part, not transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'N', 'U' );

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
      test_ = "Column-major unitriangular LSE substitution (multiple right-hand sides, lower part, transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'T', 'U' );

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
      test_ = "Column-major unitriangular LSE substitution (multiple right-hand sides, lower part, conjugate transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'C', 'U' );

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
      test_ = "Column-major triangular LSE substitution (multiple right-hand sides, upper part, not transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'N', 'N' );

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
      test_ = "Column-major triangular LSE substitution (multiple right-hand sides, upper part, transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'T', 'N' );

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
      test_ = "Column-major triangular LSE substitution (multiple right-hand sides, upper part, conjugate transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'C', 'N' );

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
      test_ = "Column-major unitriangular LSE substitution (multiple right-hand sides, upper part, not transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'U', 'N', 'U' );

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
      test_ = "Column-major unitriangular LSE substitution (multiple right-hand sides, upper part, transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'T', 'U' );

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
      test_ = "Column-major unitriangular LSE substitution (multiple right-hand sides, upper part, conjugate transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticMatrix<Type,3UL,6UL,blaze::columnMajor> B, X;
      randomize( B );

      X = B;

      blaze::trtrs( LU, X, 'L', 'C', 'U' );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the general linear system solver functions (gesv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the general linear system solver functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testGesv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major general LSE (single right-hand side)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::gesv( LU, x, ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general LSE (multiple right-hand sides)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::gesv( LU, X, ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major general LSE (single right-hand side)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::gesv( LU, x, ipiv.data() );

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
      test_ = "Column-major general LSE (multiple right-hand sides)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::gesv( LU, X, ipiv.data() );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the symmetric indefinite linear system solver functions (sysv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the symmetric indefinite linear system solver functions for
// various data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testSysv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major symmetric indefinite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::sysv( LU, x, 'L', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::sysv( LU, x, 'U', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::sysv( LU, X, 'L', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::sysv( LU, X, 'U', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major symmetric indefinite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::sysv( LU, x, 'L', ipiv.data() );

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
      test_ = "Column-major symmetric indefinite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::sysv( LU, x, 'U', ipiv.data() );

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
      test_ = "Column-major symmetric indefinite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::sysv( LU, X, 'L', ipiv.data() );

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
      test_ = "Column-major symmetric indefinite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::sysv( LU, X, 'U', ipiv.data() );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Hermitian indefinite linear system solver functions (hesv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Hermitian indefinite linear system solver functions for
// various data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testHesv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Hermitian indefinite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::hesv( LU, x, 'L', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::hesv( LU, x, 'U', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::hesv( LU, X, 'L', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::hesv( LU, X, 'U', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Hermitian indefinite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::hesv( LU, x, 'L', ipiv.data() );

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
      test_ = "Column-major Hermitian indefinite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::hesv( LU, x, 'U', ipiv.data() );

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
      test_ = "Column-major Hermitian indefinite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::hesv( LU, X, 'L', ipiv.data() );

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
      test_ = "Column-major Hermitian indefinite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<int,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::hesv( LU, X, 'U', ipiv.data() );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the positive definite linear system solver functions (posv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the positive definite linear system solver functions for
// various data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testPosv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major positive definite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::posv( LU, x, 'L' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::posv( LU, x, 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::posv( LU, X, 'L' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::posv( LU, X, 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major positive definite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::posv( LU, x, 'L' );

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
      test_ = "Column-major positive definite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::posv( LU, x, 'U' );

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
      test_ = "Column-major positive definite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::posv( LU, X, 'L' );

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
      test_ = "Column-major positive definite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::posv( LU, X, 'U' );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the triangular linear system solver functions (trsv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the triangular linear system solver functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testTrsv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major triangular LSE (single right-hand side, lower part, not transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'N', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, lower part, transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'T', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, lower part, conjugate transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'C', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, lower part, not transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'N', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, lower part, transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'T', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, lower part, conjugate transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'C', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, upper part, not transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'N', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, upper part, transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'T', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, upper part, conjugate transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'C', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, upper part, not transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'N', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, upper part, transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'T', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, upper part, conjugate transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'C', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major triangular LSE (single right-hand side, lower part, not transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'N', 'N' );

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
      test_ = "Column-major triangular LSE (single right-hand side, lower part, transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'T', 'N' );

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
      test_ = "Column-major triangular LSE (single right-hand side, lower part, conjugate transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'C', 'N' );

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
      test_ = "Column-major unitriangular LSE (single right-hand side, lower part, not transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'N', 'U' );

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
      test_ = "Column-major unitriangular LSE (single right-hand side, lower part, transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'T', 'U' );

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
      test_ = "Column-major unitriangular LSE (single right-hand side, lower part, conjugate transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'C', 'U' );

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
      test_ = "Column-major triangular LSE (single right-hand side, upper part, not transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'N', 'N' );

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
      test_ = "Column-major triangular LSE (single right-hand side, upper part, transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'T', 'N' );

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
      test_ = "Column-major triangular LSE (single right-hand side, upper part, conjugate transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'C', 'N' );

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
      test_ = "Column-major unitriangular LSE (single right-hand side, upper part, not transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'N', 'U' );

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
      test_ = "Column-major unitriangular LSE (single right-hand side, upper part, transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'T', 'U' );

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
      test_ = "Column-major unitriangular LSE (single right-hand side, upper part, conjugate transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'C', 'U' );

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

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the QR decomposition functions (geqrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the QR decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testGeqrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "QR decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      if( A != trans( B ) || tauA != tauB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      if( A != trans( B ) || tauA != tauB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QR decomposition (orgqr).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QR decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testOrgqr()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QR decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::orgqr( A, tauA.data() );
      blaze::orgqr( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::orgqr( A, tauA.data() );
      blaze::orgqr( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QR decomposition (ungqr).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QR decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testUngqr()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QR decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ungqr( A, tauA.data() );
      blaze::ungqr( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ungqr( A, tauA.data() );
      blaze::ungqr( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the RQ decomposition functions (gerqf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the RQ decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testGerqf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "RQ decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      if( A != trans( B ) || tauA != tauB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: RQ decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      if( A != trans( B ) || tauA != tauB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: RQ decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a RQ decomposition (orgrq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a RQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testOrgrq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a RQ decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::orgrq( A, tauA.data() );
      blaze::orgrq( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::orgrq( A, tauA.data() );
      blaze::orgrq( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a RQ decomposition (ungrq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a RQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testUngrq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a RQ decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ungrq( A, tauA.data() );
      blaze::ungrq( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ungrq( A, tauA.data() );
      blaze::ungrq( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the QL decomposition functions (geqlf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the QL decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testGeqlf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "QL decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      if( A != trans( B ) || tauA != tauB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QL decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      if( A != trans( B ) || tauA != tauB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QL decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QL decomposition (orgql).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QL decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testOrgql()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QL decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::orgql( A, tauA.data() );
      blaze::orgql( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::orgql( A, tauA.data() );
      blaze::orgql( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QL decomposition (ungql).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QL decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testUngql()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QL decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ungql( A, tauA.data() );
      blaze::ungql( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ungql( A, tauA.data() );
      blaze::ungql( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functions (gelqf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the LQ decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testGelqf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "LQ decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      if( A != trans( B ) || tauA != tauB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LQ decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      if( A != trans( B ) || tauA != tauB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LQ decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a LQ decomposition (orglq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a LQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testOrglq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a LQ decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::orglq( A, tauA.data() );
      blaze::orglq( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::orglq( A, tauA.data() );
      blaze::orglq( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a LQ decomposition (unglq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a LQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void OperationTest::testUnglq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a LQ decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::unglq( A, tauA.data() );
      blaze::unglq( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::unglq( A, tauA.data() );
      blaze::unglq( B, tauB.data() );

      if( A != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
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
/*!\brief Testing the LAPACK functionality.
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
/*!\brief Macro for the execution of the LAPACK operation test.
*/
#define RUN_LAPACK_OPERATION_TEST \
   blazetest::mathtest::lapack::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lapack

} // namespace mathtest

} // namespace blazetest

#endif
