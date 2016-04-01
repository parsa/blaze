//=================================================================================================
/*!
//  \file blazetest/mathtest/inversion/DenseTest.h
//  \brief Header file for the dense matrix inversion test
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

#ifndef _BLAZETEST_MATHTEST_INVERSION_DENSETEST_H_
#define _BLAZETEST_MATHTEST_INVERSION_DENSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Aliases.h>
#include <blaze/math/DenseMatrix.h>
#include <blaze/math/DenseSubmatrix.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/util/Random.h>
#include <blazetest/system/LAPACK.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace inversion {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all dense matrix inversion tests.
//
// This class represents a test suite for the dense matrix inversion functionality. It performs
// a series of matrix inversions on all dense matrix types of the Blaze library.
*/
class DenseTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit DenseTest();
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
   void testSpecific();

   template< typename Type >
   void testRandom1x1();

   template< typename Type >
   void testRandom2x2();

   template< typename Type >
   void testRandom3x3();

   template< typename Type >
   void testRandom4x4();

   template< typename Type >
   void testRandom5x5();

   template< typename Type >
   void testRandom6x6();

   template< typename Type >
   void testRandomNxN();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename MT, bool SO >
   void initializeForLU( blaze::DenseMatrix<MT,SO>& matrix );

   template< typename MT, bool SO >
   void initializeForLDLT( blaze::DenseMatrix<MT,SO>& matrix );

   template< typename MT, bool SO >
   void initializeForLDLH( blaze::DenseMatrix<MT,SO>& matrix );

   template< typename MT, bool SO >
   void initializeForLLH( blaze::DenseMatrix<MT,SO>& matrix );
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
/*!\brief Test of the inversion functionality with random 1x1 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix inversion for random 1x1 matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom1x1()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Matrix inversion (1x1)";

   Type A;
   resize( A, 1UL, 1UL );

   do {
      randomize( A );
   } while( blaze::isDefault( det( A ) ) );

   Type B( A );

   blaze::invert( A );
   blaze::invertNxN<blaze::byLU>( B );

   if( A != B ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Matrix inversion failed\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( blaze::ElementType_<Type> ).name() << "\n"
          << "   Result invert   ():\n" << A << "\n"
          << "   Result invertNxN():\n" << B << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the inversion functionality with random 2x2 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix inversion for random 2x2 matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom2x2()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Matrix inversion (2x2)";

   Type A;
   resize( A, 2UL, 2UL );

   do {
      randomize( A );
   } while( blaze::isDefault( det( A ) ) );

   Type B( A );
   Type C( A );

   blaze::invert( A );
   blaze::invert2x2( B );
   blaze::invertNxN<blaze::byLU>( C );

   if( A != B || A != C ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Matrix inversion failed\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( blaze::ElementType_<Type> ).name() << "\n"
          << "   Result invert   ():\n" << A << "\n"
          << "   Result invert2x2():\n" << B << "\n"
          << "   Result invertNxN():\n" << C << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the inversion functionality with random 3x3 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix inversion for random 3x3 matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom3x3()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Matrix inversion (3x3)";

   Type A;
   resize( A, 3UL, 3UL );

   do {
      randomize( A );
   } while( blaze::isDefault( det( A ) ) );

   Type B( A );
   Type C( A );

   blaze::invert( A );
   blaze::invert3x3( B );
   blaze::invertNxN<blaze::byLU>( C );

   if( A != B || A != C ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Matrix inversion failed\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( blaze::ElementType_<Type> ).name() << "\n"
          << "   Result invert   ():\n" << A << "\n"
          << "   Result invert3x3():\n" << B << "\n"
          << "   Result invertNxN():\n" << C << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the inversion functionality with random 4x4 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix inversion for random 4x4 matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom4x4()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Matrix inversion (4x4)";

   Type A;
   resize( A, 4UL, 4UL );

   do {
      randomize( A );
   } while( blaze::isDefault( det( A ) ) );

   Type B( A );
   Type C( A );

   blaze::invert( A );
   blaze::invert4x4( B );
   blaze::invertNxN<blaze::byLU>( C );

   if( A != B || A != C ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Matrix inversion failed\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( blaze::ElementType_<Type> ).name() << "\n"
          << "   Result invert   ():\n" << A << "\n"
          << "   Result invert4x4():\n" << B << "\n"
          << "   Result invertNxN():\n" << C << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the inversion functionality with random 5x5 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix inversion for random 5x5 matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom5x5()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Matrix inversion (5x5)";

   Type A;
   resize( A, 5UL, 5UL );

   do {
      randomize( A );
   } while( blaze::isDefault( det( A ) ) );

   Type B( A );
   Type C( A );

   blaze::invert( A );
   blaze::invert5x5( B );
   blaze::invertNxN<blaze::byLU>( C );

   if( A != B || A != C ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Matrix inversion failed\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( blaze::ElementType_<Type> ).name() << "\n"
          << "   Result invert   ():\n" << A << "\n"
          << "   Result invert5x5():\n" << B << "\n"
          << "   Result invertNxN():\n" << C << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the inversion functionality with random 6x6 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix inversion for random 6x6 matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom6x6()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Matrix inversion (6x6)";

   Type A;
   resize( A, 6UL, 6UL );

   do {
      randomize( A );
   } while( blaze::isDefault( det( A ) ) );

   Type B( A );
   Type C( A );

   blaze::invert( A );
   blaze::invert6x6( B );
   blaze::invertNxN<blaze::byLU>( C );

   if( A != B || A != C ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Matrix inversion failed\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( blaze::ElementType_<Type> ).name() << "\n"
          << "   Result invert   ():\n" << A << "\n"
          << "   Result invert6x6():\n" << B << "\n"
          << "   Result invertNxN():\n" << C << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the inversion functionality with random \f$ N \times N \f$ matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix inversion for random \f$ N \times N \f$ matrices. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandomNxN()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::invert;
   using blaze::byDefault;
   using blaze::byLU;
   using blaze::byLDLT;
   using blaze::byLDLH;
   using blaze::byLLH;

   typedef blaze::ElementType_<Type>  ET;


   //=====================================================================================
   // Matrix default inversion
   //=====================================================================================

   {
      test_ = "Matrix inversion (default)";

      Type A;
      initializeForLU( A );
      Type B( A );

      invert<byDefault>( B );

      if( !isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Submatrix inversion (default)";

      Type A;
      initializeForLU( A );
      Type B( A );

      blaze::DenseSubmatrix<Type> sub( B, 0UL, 0UL, A.rows(), A.columns() );
      invert<byDefault>( sub );

      if( !isIdentity( A * sub ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix inversion by LU decomposition
   //=====================================================================================

   {
      test_ = "Matrix inversion (LU)";

      Type A;
      initializeForLU( A );
      Type B( A );

      invert<byLU>( B );

      if( !isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Submatrix inversion (LU)";

      Type A;
      initializeForLU( A );
      Type B( A );

      blaze::DenseSubmatrix<Type> sub( B, 0UL, 0UL, A.rows(), A.columns() );
      invert<byLU>( sub );

      if( !isIdentity( A * sub ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix inversion by LDLT (Bunch-Kaufman) decomposition
   //=====================================================================================

   {
      test_ = "Matrix inversion (LDLT/Bunch-Kaufman)";

      Type A;
      initializeForLDLT( A );
      Type B( A );

      invert<byLDLT>( B );

      if( !isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Submatrix inversion (LDLT/Bunch-Kaufman)";

      Type A;
      initializeForLDLT( A );
      Type B( A );

      blaze::DenseSubmatrix<Type> sub( B, 0UL, 0UL, A.rows(), A.columns() );
      invert<byLDLT>( sub );

      if( !isIdentity( A * sub ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix inversion by LDLH (Bunch-Kaufman) decomposition
   //=====================================================================================

   {
      test_ = "Matrix inversion (LDLH/Bunch-Kaufman)";

      Type A;
      initializeForLDLH( A );
      Type B( A );

      invert<byLDLH>( B );

      if( !isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Submatrix inversion (LDLH/Bunch-Kaufman)";

      Type A;
      initializeForLDLH( A );
      Type B( A );

      blaze::DenseSubmatrix<Type> sub( B, 0UL, 0UL, A.rows(), A.columns() );
      invert<byLDLH>( sub );

      if( !isIdentity( A * sub ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix inversion by LLH (Cholesky) decomposition
   //=====================================================================================

   {
      test_ = "Matrix inversion (LLH/Cholesky)";

      Type A;
      initializeForLLH( A );
      Type B( A );

      invert<byLLH>( B );

      if( !isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Submatrix inversion (LLH/Cholesky)";

      Type A;
      initializeForLLH( A );
      Type B( A );

      blaze::DenseSubmatrix<Type> sub( B, 0UL, 0UL, A.rows(), A.columns() );
      invert<byLLH>( sub );

      if( !isIdentity( A * sub ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix (A):\n" << A << "\n"
             << "   Result (B):\n" << B << "\n"
             << "   A * B =\n" << ( A * B ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of the given dense matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT, bool SO >
void DenseTest::initializeForLU( blaze::DenseMatrix<MT,SO>& matrix )
{
   const size_t size( blaze::rand<size_t>( 7UL, 14UL ) );

   resize( ~matrix, size, size );
   randomize( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given dense matrix for a LDLT-based matrix inversion.
//
// \return void
*/
template< typename MT, bool SO >
void DenseTest::initializeForLDLT( blaze::DenseMatrix<MT,SO>& matrix )
{
   const size_t size( blaze::rand<size_t>( 7UL, 14UL ) );

   resize( ~matrix, size, size );
   makeSymmetric( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given dense matrix for a LDLH-based matrix inversion.
//
// \return void
*/
template< typename MT, bool SO >
void DenseTest::initializeForLDLH( blaze::DenseMatrix<MT,SO>& matrix )
{
   const size_t size( blaze::rand<size_t>( 7UL, 14UL ) );

   resize( ~matrix, size, size );
   makeHermitian( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given dense matrix for a LLH-based matrix inversion.
//
// \return void
*/
template< typename MT, bool SO >
void DenseTest::initializeForLLH( blaze::DenseMatrix<MT,SO>& matrix )
{
   const size_t size( blaze::rand<size_t>( 7UL, 14UL ) );

   resize( ~matrix, size, size );
   makePositiveDefinite( ~matrix );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the dense matrix inversion.
//
// \return void
*/
void runTest()
{
   DenseTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the dense matrix inversion test.
*/
#define RUN_INVERSION_DENSE_TEST \
   blazetest::mathtest::inversion::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace inversion

} // namespace mathtest

} // namespace blazetest

#endif
