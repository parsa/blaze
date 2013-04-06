//=================================================================================================
/*!
//  \file src/mathtest/smatsmatmult/MCbMCa.cpp
//  \brief Source file for the MCbMCa sparse matrix/sparse matrix multiplication math test
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
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/smatsmatmult/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'MCbMCa'..." << std::endl;

   using blazetest::mathtest::TypeA;
   using blazetest::mathtest::TypeB;

   try
   {
      // Matrix type definitions
      typedef blaze::CompressedMatrix<TypeB>  MCb;
      typedef blaze::CompressedMatrix<TypeA>  MCa;

      // Creator type definitions
      typedef blazetest::Creator<MCb>  CMCb;
      typedef blazetest::Creator<MCa>  CMCa;

      // Running tests with small matrices
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=6UL; ++j ) {
            for( size_t k=0UL; k<=6UL; ++k ) {
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j,     0UL ), CMCa( j, k,     0UL ) );
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j,     0UL ), CMCa( j, k, 0.3*j*k ) );
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j,     0UL ), CMCa( j, k,     j*k ) );
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j, 0.3*i*j ), CMCa( j, k,     0UL ) );
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j, 0.3*i*j ), CMCa( j, k, 0.3*j*k ) );
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j, 0.3*i*j ), CMCa( j, k,     j*k ) );
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j,     i*j ), CMCa( j, k,     0UL ) );
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j,     i*j ), CMCa( j, k, 0.3*j*k ) );
               RUN_SMATSMATMULT_OPERATION_TEST( CMCb( i, j,     i*j ), CMCa( j, k,     j*k ) );
            }
         }
      }

      // Running tests with large matrices
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 15UL,  37UL,  7UL ), CMCa(  37UL, 15UL,  7UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 15UL,  37UL,  7UL ), CMCa(  37UL, 63UL, 13UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 37UL,  37UL,  7UL ), CMCa(  37UL, 37UL,  7UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 63UL,  37UL, 13UL ), CMCa(  37UL, 15UL,  7UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 63UL,  37UL, 13UL ), CMCa(  37UL, 63UL, 13UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 16UL,  32UL,  8UL ), CMCa(  32UL, 16UL,  8UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 16UL,  32UL,  8UL ), CMCa(  32UL, 64UL, 16UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 32UL,  32UL,  8UL ), CMCa(  32UL, 32UL,  8UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 64UL,  32UL, 16UL ), CMCa(  32UL, 16UL,  8UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCb( 64UL,  32UL, 16UL ), CMCa(  32UL, 64UL, 16UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse matrix/sparse matrix multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
