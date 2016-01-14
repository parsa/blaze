//=================================================================================================
/*!
//  \file blazetest/mathtest/lu/DenseTest.h
//  \brief Header file for the dense matrix LU test
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

#ifndef _BLAZETEST_MATHTEST_LU_DENSETEST_H_
#define _BLAZETEST_MATHTEST_LU_DENSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blazetest/system/LAPACK.h>


namespace blazetest {

namespace mathtest {

namespace lu {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all dense matrix LU tests.
//
// This class represents a test suite for the dense matrix LU decomposition functionality. It
// performs a series of LU decompositions on all dense matrix types of the Blaze library.
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
   template< typename Type > void testGeneral();
   template< typename Type > void testSymmetric();
   template< typename Type > void testHermitian();
   template< typename Type > void testLower();
   template< typename Type > void testUniLower();
   template< typename Type > void testUpper();
   template< typename Type > void testUniUpper();
   template< typename Type > void testDiagonal();
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
/*!\brief Test of the LU decomposition functionality for general matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for general matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testGeneral()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major general matrix (3x3)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general matrix (2x5)";

      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,2UL,2UL,blaze::rowMajor> > L;
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> U;
      blaze::StaticMatrix<Type,5UL,5UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general matrix (5x2)";

      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,2UL,2UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,2UL,2UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major general matrix (3x3)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major general matrix (2x5)";

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,2UL,2UL,blaze::columnMajor> > L;
      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> U;
      blaze::StaticMatrix<Type,2UL,2UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major general matrix (5x2)";

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,2UL,2UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,5UL,5UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LU decomposition functionality for symmetric matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for symmetric matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testSymmetric()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major symmetric matrix (3x3)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major symmetric matrix (3x3)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LU decomposition functionality for Hermitian matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for Hermitian matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testHermitian()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Hermitian matrix (3x3)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Hermitian matrix (3x3)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LU decomposition functionality for lower matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for lower matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testLower()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major lower matrix (3x3)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major lower matrix (3x3)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LU decomposition functionality for unilower matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for unilower matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testUniLower()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major unilower matrix (3x3)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major unilower matrix (3x3)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LU decomposition functionality for upper matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for upper matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testUpper()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major upper matrix (3x3)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major upper matrix (3x3)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LU decomposition functionality for uniupper matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for uniupper matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testUniUpper()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major uniupper matrix (3x3)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major uniupper matrix (3x3)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LU decomposition functionality for diagonal matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for diagonal matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testDiagonal()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major diagonal matrix (3x3)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LUP( L*U*P );

      if( LUP != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << LUP << "\n"
             << "   Expected result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major diagonal matrix (3x3)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > L;
      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > U;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> P;

      blaze::lu( A, L, U, P );

      const blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> PLU( P*L*U );

      if( PLU != A ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Result:\n" << PLU << "\n"
             << "   Expected result:\n" << A << "\n";
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
/*!\brief Testing the dense matrix LU decomposition.
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
/*!\brief Macro for the execution of the dense matrix LU test.
*/
#define RUN_LU_DENSE_TEST \
   blazetest::mathtest::lu::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lu

} // namespace mathtest

} // namespace blazetest

#endif
