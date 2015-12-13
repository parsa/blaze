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
#include <blaze/math/DenseSubmatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blazetest/system/LAPACK.h>


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
   template< typename Type >
   void testInversion();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename MT, bool SO > void initializeForLU( blaze::DenseMatrix<MT,SO>& matrix );
   template< typename MT >          void initializeForLU( blaze::SymmetricMatrix<MT>& sym );
   template< typename MT >          void initializeForLU( blaze::HermitianMatrix<MT>& herm );
   template< typename MT >          void initializeForLU( blaze::LowerMatrix<MT>& lower );
   template< typename MT >          void initializeForLU( blaze::UniLowerMatrix<MT>& lower );
   template< typename MT >          void initializeForLU( blaze::UpperMatrix<MT>& upper );
   template< typename MT >          void initializeForLU( blaze::UniUpperMatrix<MT>& upper );
   template< typename MT >          void initializeForLU( blaze::DiagonalMatrix<MT>& diag );

   template< typename MT, bool SO > void initializeForCholesky( blaze::DenseMatrix<MT,SO>& matrix );
   template< typename MT >          void initializeForCholesky( blaze::SymmetricMatrix<MT>& sym );
   template< typename MT >          void initializeForCholesky( blaze::HermitianMatrix<MT>& herm );
   template< typename MT >          void initializeForCholesky( blaze::LowerMatrix<MT>& lower );
   template< typename MT >          void initializeForCholesky( blaze::UniLowerMatrix<MT>& lower );
   template< typename MT >          void initializeForCholesky( blaze::UpperMatrix<MT>& upper );
   template< typename MT >          void initializeForCholesky( blaze::UniUpperMatrix<MT>& upper );
   template< typename MT >          void initializeForCholesky( blaze::DiagonalMatrix<MT>& diag );
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
/*!\brief Test of the QR decomposition functionality.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the invert() function for a specific matrix type. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testInversion()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::invert;
   using blaze::byLU;
   using blaze::byCholesky;
   using blaze::equal;

   typedef typename Type::ElementType  ET;


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
             << "   Initial matrix:\n" << A << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Submatrix inversion (LU)";

      Type A;
      initializeForLU( A );
      Type B( A );

      blaze::DenseSubmatrix<Type> sub( B, 0UL, 0UL, 3UL, 3UL );
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
             << "   Initial matrix:\n" << A << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix inversion by Cholesky decomposition
   //=====================================================================================

   {
      test_ = "Matrix inversion (Cholesky)";

      Type A;
      initializeForCholesky( A );
      Type B( A );

      invert<byCholesky>( B );

      if( !isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix:\n" << A << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Submatrix inversion (Cholesky)";

      Type A;
      initializeForCholesky( A );
      Type B( A );

      blaze::DenseSubmatrix<Type> sub( B, 0UL, 0UL, 3UL, 3UL );
      invert<byCholesky>( sub );

      if( !isIdentity( A * sub ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix inversion failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ET ).name() << "\n"
             << "   Initial matrix:\n" << A << "\n"
             << "   Result:\n" << B << "\n";
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
/*!\brief Initialization of the given matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT, bool SO >
void DenseTest::initializeForLU( blaze::DenseMatrix<MT,SO>& matrix )
{
   typedef typename MT::ElementType  ET;

   resize( ~matrix, 3UL, 3UL );
   reset( ~matrix );

   (~matrix)(0,0) = ET(1);
   (~matrix)(1,1) = ET(1);
   (~matrix)(2,0) = ET(1);
   (~matrix)(2,1) = ET(1);
   (~matrix)(2,2) = ET(1);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given symmetric matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForLU( blaze::SymmetricMatrix<MT>& sym )
{
   typedef typename blaze::SymmetricMatrix<MT>::ElementType  ET;

   resize( sym, 3UL, 3UL );
   reset( sym );

   sym(0,0) = ET(1);
   sym(0,2) = ET(1);
   sym(1,1) = ET(1);
   sym(1,2) = ET(1);
   sym(2,2) = ET(1);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given Hermitian matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForLU( blaze::HermitianMatrix<MT>& herm )
{
   typedef typename blaze::HermitianMatrix<MT>::ElementType  ET;

   resize( herm, 3UL, 3UL );
   reset( herm );

   herm(0,0) = ET(1);
   herm(0,2) = ET(1);
   herm(1,1) = ET(1);
   herm(1,2) = ET(1);
   herm(2,2) = ET(1);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given lower matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForLU( blaze::LowerMatrix<MT>& lower )
{
   typedef typename blaze::LowerMatrix<MT>::ElementType  ET;

   resize( lower, 3UL, 3UL );
   reset( lower );

   lower(0,0) = ET(1);
   lower(1,1) = ET(1);
   lower(2,0) = ET(1);
   lower(2,1) = ET(1);
   lower(2,2) = ET(1);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given unilower matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForLU( blaze::UniLowerMatrix<MT>& lower )
{
   typedef typename blaze::UniLowerMatrix<MT>::ElementType  ET;

   resize( lower, 3UL, 3UL );
   reset( lower );

   lower(2,0) = ET(1);
   lower(2,1) = ET(1);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given upper matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForLU( blaze::UpperMatrix<MT>& upper )
{
   typedef typename blaze::UpperMatrix<MT>::ElementType  ET;

   resize( upper, 3UL, 3UL );
   reset( upper );

   upper(0,0) = ET(1);
   upper(0,2) = ET(1);
   upper(1,1) = ET(1);
   upper(1,2) = ET(1);
   upper(2,2) = ET(1);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given uniupper matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForLU( blaze::UniUpperMatrix<MT>& upper )
{
   typedef typename blaze::UniUpperMatrix<MT>::ElementType  ET;

   resize( upper, 3UL, 3UL );
   reset( upper );

   upper(0,2) = ET(1);
   upper(1,2) = ET(1);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given diagonal matrix for a LU-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForLU( blaze::DiagonalMatrix<MT>& diag )
{
   typedef typename blaze::DiagonalMatrix<MT>::ElementType  ET;

   resize( diag, 3UL, 3UL );
   reset( diag );

   diag(0,0) = ET(4);
   diag(1,1) = ET(4);
   diag(2,2) = ET(4);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given matrix for a Cholesky-based matrix inversion.
//
// \return void
*/
template< typename MT, bool SO >
void DenseTest::initializeForCholesky( blaze::DenseMatrix<MT,SO>& matrix )
{
   typedef typename MT::ElementType  ET;

   resize( ~matrix, 3UL, 3UL );

   (~matrix)(0,0) = ET(1);
   (~matrix)(0,1) = ET(1);
   (~matrix)(0,2) = ET(1);
   (~matrix)(1,0) = ET(1);
   (~matrix)(1,1) = ET(2);
   (~matrix)(1,2) = ET(2);
   (~matrix)(2,0) = ET(1);
   (~matrix)(2,1) = ET(2);
   (~matrix)(2,2) = ET(4);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given symmetric matrix for a Cholesky-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForCholesky( blaze::SymmetricMatrix<MT>& sym )
{
   typedef typename blaze::SymmetricMatrix<MT>::ElementType  ET;

   resize( sym, 3UL, 3UL );
   reset( sym );

   sym(0,0) = ET(1);
   sym(0,1) = ET(1);
   sym(0,2) = ET(1);
   sym(1,1) = ET(2);
   sym(1,2) = ET(2);
   sym(2,2) = ET(4);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given Hermitian matrix for a Cholesky-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForCholesky( blaze::HermitianMatrix<MT>& herm )
{
   typedef typename blaze::HermitianMatrix<MT>::ElementType  ET;

   resize( herm, 3UL, 3UL );
   reset( herm );

   herm(0,0) = ET(1);
   herm(0,1) = ET(1);
   herm(0,2) = ET(1);
   herm(1,1) = ET(2);
   herm(1,2) = ET(2);
   herm(2,2) = ET(4);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given lower matrix for a Cholesky-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForCholesky( blaze::LowerMatrix<MT>& lower )
{
   typedef typename blaze::LowerMatrix<MT>::ElementType  ET;

   resize( lower, 3UL, 3UL );
   reset( lower );

   lower(0,0) = ET(4);
   lower(1,1) = ET(4);
   lower(2,2) = ET(4);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given unilower matrix for a Cholesky-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForCholesky( blaze::UniLowerMatrix<MT>& lower )
{
   typedef typename blaze::UniLowerMatrix<MT>::ElementType  ET;

   resize( lower, 3UL, 3UL );
   reset( lower );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given upper matrix for a Cholesky-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForCholesky( blaze::UpperMatrix<MT>& upper )
{
   typedef typename blaze::UpperMatrix<MT>::ElementType  ET;

   resize( upper, 3UL, 3UL );
   reset( upper );

   upper(0,0) = ET(4);
   upper(1,1) = ET(4);
   upper(2,2) = ET(4);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given uniupper matrix for a Cholesky-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForCholesky( blaze::UniUpperMatrix<MT>& upper )
{
   typedef typename blaze::UniUpperMatrix<MT>::ElementType  ET;

   resize( upper, 3UL, 3UL );
   reset( upper );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of the given diagonal matrix for a Cholesky-based matrix inversion.
//
// \return void
*/
template< typename MT >
void DenseTest::initializeForCholesky( blaze::DiagonalMatrix<MT>& diag )
{
   typedef typename blaze::DiagonalMatrix<MT>::ElementType  ET;

   resize( diag, 3UL, 3UL );
   reset( diag );

   diag(0,0) = ET(4);
   diag(1,1) = ET(4);
   diag(2,2) = ET(4);
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
