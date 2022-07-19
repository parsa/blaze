//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/svd/DenseTest.h
//  \brief Header file for the dense matrix singular value test
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_SVD_DENSETEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_SVD_DENSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <string>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace svd {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all dense matrix singular value tests.
//
// This class represents a test suite for the dense matrix singular value functionality. It
// performs a series of singular value computations on several dense matrix types of the Blaze
// library.
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
   template< typename Type, bool SO >
   void testMatrixRandom();
   template< typename Type, bool SO >
   void testMatrixRandomSingle(size_t m, size_t n, bool square);

   void testGeneral();
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
/*!\brief Test of the SVD decomposition with a randomly initialized matrix of the given type.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix SVD decomposition for a randomly initialized matrix of the
// given type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type, bool SO >
void DenseTest::testMatrixRandom()
{
   size_t m( blaze::rand<size_t>( 4UL, 8UL ) );
   size_t n( m );

   // test m == n
   testMatrixRandomSingle<Type, SO>(m, n, false);
   
   if (!blaze::IsSquare<Type>::value) {
      // test m > n, squareV = false
      n = blaze::rand<size_t>( 2UL, m - 1 );
      testMatrixRandomSingle<Type, SO>(m, n, false);

      // test m > n, squareV = true
      testMatrixRandomSingle<Type, SO>(m, n, true);

      // test m < n, squareV = false
      n = blaze::rand<size_t>( m + 1, 10UL );
      testMatrixRandomSingle<Type, SO>(m, n, false);

      // test m < n, squareV = true
      testMatrixRandomSingle<Type, SO>(m, n, true);
   }
}

//*************************************************************************************************
/*!\brief Test of the QR decomposition with a randomly initialized matrix of the given type.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix QR decomposition for a randomly initialized matrix of the
// given type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type, bool SO >
void DenseTest::testMatrixRandomSingle(size_t m, size_t n, bool square)
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using MT = blaze::RemoveAdaptor_t<Type>;

   DynamicMatrix<Type, SO> A( m, n );
   randomize( A );

   DynamicMatrix<Type, SO> U, V;
   DynamicVector<Type, rowVector> s;

   svd(A, U, s, V, square);

   DynamicMatrix<Type, SO> S( U.columns(), V.rows(), 0 );
   diagonal(S) = s;

   DynamicMatrix<Type, SO> USV = U * S * V;

   if (square) {
      if (U.rows() != U.columns() || U.rows() != m) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
               << " Error: Singular value computation failed\n"
               << " Details:\n"
               << "   Random seed = " << blaze::getSeed() << "\n"
               << "   U # of rows:\n" << U.rows() << "\n"
               << "   U # of columns:\n" << U.columns() << "\n"
               << "   Expected # of rows/columns:\n" << m << "\n";
         throw std::runtime_error( oss.str() );
      }
      if (V.rows() != V.columns() || V.rows() != n) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
               << " Error: Singular value computation failed\n"
               << " Details:\n"
               << "   Random seed = " << blaze::getSeed() << "\n"
               << "   V # of rows:\n" << V.rows() << "\n"
               << "   V # of columns:\n" << V.columns() << "\n"
               << "   Expected # of rows/columns:\n" << n << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   if (A != USV) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
            << " Error: Singular value computation failed\n"
            << " Details:\n"
            << "   Random seed = " << blaze::getSeed() << "\n"
            << "   singular values:\n" << s << "\n"
            << "   left singular vectors:\n" << U << "\n"
            << "   right singular vectors:\n" << V << "\n"
            << "   Product:\n" << USV << "\n"
            << "   Expected Result:\n" << A << "\n";
      throw std::runtime_error( oss.str() );
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
/*!\brief Testing the dense matrix singular value functionality.
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
/*!\brief Macro for the execution of the dense matrix singular value test.
*/
#define RUN_DENSE_SVD_TEST \
   blazetest::mathtest::operations::svd::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace svd

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
