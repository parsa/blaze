//=================================================================================================
/*!
//  \file src/mathtest/inversion/DenseTest.cpp
//  \brief Source file for the dense matrix inversion test
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
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/inversion/DenseTest.h>


namespace blazetest {

namespace mathtest {

namespace inversion {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest inversion test.
//
// \exception std::runtime_error Matrix inversion error detected.
*/
DenseTest::DenseTest()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::HermitianMatrix;
   using blaze::LowerMatrix;
   using blaze::UniLowerMatrix;
   using blaze::UpperMatrix;
   using blaze::UniUpperMatrix;
   using blaze::DiagonalMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   typedef blaze::complex<float>   cfloat;
   typedef blaze::complex<double>  cdouble;


   //=====================================================================================
   // Specific matrix tests
   //=====================================================================================

   testSpecific();


   //=====================================================================================
   // Random 1x1 matrix tests
   //=====================================================================================

   testRandom1x1< DynamicMatrix<double ,rowMajor> >();
   testRandom1x1< DynamicMatrix<cdouble,rowMajor> >();

   testRandom1x1< DynamicMatrix<double ,columnMajor> >();
   testRandom1x1< DynamicMatrix<cdouble,columnMajor> >();

   testRandom1x1< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom1x1< HermitianMatrix< DynamicMatrix<cdouble,rowMajor> > >();
   testRandom1x1< LowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom1x1< UniLowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom1x1< UpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom1x1< UniUpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom1x1< DiagonalMatrix< DynamicMatrix<double,rowMajor> > >();

   testRandom1x1< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom1x1< HermitianMatrix< DynamicMatrix<cdouble,columnMajor> > >();
   testRandom1x1< LowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom1x1< UniLowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom1x1< UpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom1x1< UniUpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom1x1< DiagonalMatrix< DynamicMatrix<double,columnMajor> > >();


   //=====================================================================================
   // Random 2x2 matrix tests
   //=====================================================================================

   testRandom2x2< DynamicMatrix<double ,rowMajor> >();
   testRandom2x2< DynamicMatrix<cdouble,rowMajor> >();

   testRandom2x2< DynamicMatrix<double ,columnMajor> >();
   testRandom2x2< DynamicMatrix<cdouble,columnMajor> >();

   testRandom2x2< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom2x2< HermitianMatrix< DynamicMatrix<cdouble,rowMajor> > >();
   testRandom2x2< LowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom2x2< UniLowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom2x2< UpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom2x2< UniUpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom2x2< DiagonalMatrix< DynamicMatrix<double,rowMajor> > >();

   testRandom2x2< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom2x2< HermitianMatrix< DynamicMatrix<cdouble,columnMajor> > >();
   testRandom2x2< LowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom2x2< UniLowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom2x2< UpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom2x2< UniUpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom2x2< DiagonalMatrix< DynamicMatrix<double,columnMajor> > >();


   //=====================================================================================
   // Random 3x3 matrix tests
   //=====================================================================================

   testRandom3x3< DynamicMatrix<double ,rowMajor> >();
   testRandom3x3< DynamicMatrix<cdouble,rowMajor> >();

   testRandom3x3< DynamicMatrix<double ,columnMajor> >();
   testRandom3x3< DynamicMatrix<cdouble,columnMajor> >();

   testRandom3x3< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom3x3< HermitianMatrix< DynamicMatrix<cdouble,rowMajor> > >();
   testRandom3x3< LowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom3x3< UniLowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom3x3< UpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom3x3< UniUpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom3x3< DiagonalMatrix< DynamicMatrix<double,rowMajor> > >();

   testRandom3x3< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom3x3< HermitianMatrix< DynamicMatrix<cdouble,columnMajor> > >();
   testRandom3x3< LowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom3x3< UniLowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom3x3< UpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom3x3< UniUpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom3x3< DiagonalMatrix< DynamicMatrix<double,columnMajor> > >();


   //=====================================================================================
   // Random 4x4 matrix tests
   //=====================================================================================

   testRandom4x4< DynamicMatrix<double ,rowMajor> >();
   testRandom4x4< DynamicMatrix<cdouble,rowMajor> >();

   testRandom4x4< DynamicMatrix<double ,columnMajor> >();
   testRandom4x4< DynamicMatrix<cdouble,columnMajor> >();

   testRandom4x4< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom4x4< HermitianMatrix< DynamicMatrix<cdouble,rowMajor> > >();
   testRandom4x4< LowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom4x4< UniLowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom4x4< UpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom4x4< UniUpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom4x4< DiagonalMatrix< DynamicMatrix<double,rowMajor> > >();

   testRandom4x4< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom4x4< HermitianMatrix< DynamicMatrix<cdouble,columnMajor> > >();
   testRandom4x4< LowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom4x4< UniLowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom4x4< UpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom4x4< UniUpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom4x4< DiagonalMatrix< DynamicMatrix<double,columnMajor> > >();


   //=====================================================================================
   // Random 5x5 matrix tests
   //=====================================================================================

   testRandom5x5< DynamicMatrix<double ,rowMajor> >();
   testRandom5x5< DynamicMatrix<cdouble,rowMajor> >();

   testRandom5x5< DynamicMatrix<double ,columnMajor> >();
   testRandom5x5< DynamicMatrix<cdouble,columnMajor> >();

   testRandom5x5< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom5x5< HermitianMatrix< DynamicMatrix<cdouble,rowMajor> > >();
   testRandom5x5< LowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom5x5< UniLowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom5x5< UpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom5x5< UniUpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom5x5< DiagonalMatrix< DynamicMatrix<double,rowMajor> > >();

   testRandom5x5< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom5x5< HermitianMatrix< DynamicMatrix<cdouble,columnMajor> > >();
   testRandom5x5< LowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom5x5< UniLowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom5x5< UpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom5x5< UniUpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom5x5< DiagonalMatrix< DynamicMatrix<double,columnMajor> > >();


   //=====================================================================================
   // Random 6x6 matrix tests
   //=====================================================================================

   testRandom6x6< DynamicMatrix<double ,rowMajor> >();
   testRandom6x6< DynamicMatrix<cdouble,rowMajor> >();

   testRandom6x6< DynamicMatrix<double ,columnMajor> >();
   testRandom6x6< DynamicMatrix<cdouble,columnMajor> >();

   testRandom6x6< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom6x6< HermitianMatrix< DynamicMatrix<cdouble,rowMajor> > >();
   testRandom6x6< LowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom6x6< UniLowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom6x6< UpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom6x6< UniUpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom6x6< DiagonalMatrix< DynamicMatrix<double,rowMajor> > >();

   testRandom6x6< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom6x6< HermitianMatrix< DynamicMatrix<cdouble,columnMajor> > >();
   testRandom6x6< LowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom6x6< UniLowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom6x6< UpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom6x6< UniUpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom6x6< DiagonalMatrix< DynamicMatrix<double,columnMajor> > >();


   //=====================================================================================
   // Random NxN matrix tests
   //=====================================================================================

   testRandomNxN< DynamicMatrix<double ,rowMajor> >();
   testRandomNxN< DynamicMatrix<cdouble,rowMajor> >();

   testRandomNxN< DynamicMatrix<double ,columnMajor> >();
   testRandomNxN< DynamicMatrix<cdouble,columnMajor> >();

   testRandomNxN< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandomNxN< HermitianMatrix< DynamicMatrix<cdouble,rowMajor> > >();
   testRandomNxN< LowerMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandomNxN< UpperMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandomNxN< DiagonalMatrix< DynamicMatrix<double,rowMajor> > >();

   testRandomNxN< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandomNxN< HermitianMatrix< DynamicMatrix<cdouble,columnMajor> > >();
   testRandomNxN< LowerMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandomNxN< UpperMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandomNxN< DiagonalMatrix< DynamicMatrix<double,columnMajor> > >();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the determinant functionality with specific, predetermined matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function computes determinants for specific, predetermined matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSpecific()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      {
         test_ = "Row-major dense matrix inversion (0x0)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A;

         invert( A );

         if( A.rows() != 0UL || A.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix inversion failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix inversion (non-square)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 2UL, 3UL );

         try {
            invert( A );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inversion of a non-square matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      {
         test_ = "Column-major dense matrix inversion (0x0)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A;

         invert( A );

         if( A.rows() != 0UL || A.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix inversion failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix inversion (non-square)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 2UL, 3UL );

         try {
            invert( A );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inversion of a non-square matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

#endif
}
//*************************************************************************************************

} // namespace inversion

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
   std::cout << "   Running dense matrix inversion test..." << std::endl;

   try
   {
      RUN_INVERSION_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix inversion test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
