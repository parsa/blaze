//=================================================================================================
/*!
//  \file src/mathtest/smatsmatadd/MCbMCa.cpp
//  \brief Source file for the MCbMCa sparse matrix/sparse matrix addition math test
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
#include <blazetest/mathtest/SMatSMatAdd.h>
#include <blazetest/system/MathTest.h>
#include <blazetest/util/Creator.h>


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
      typedef blaze::CompressedMatrix<TypeA>  MCa;
      typedef blaze::CompressedMatrix<TypeB>  MCb;

      // Creator type definitions
      typedef blazetest::Creator<MCa>  CMCa;
      typedef blazetest::Creator<MCb>  CMCb;

      // Running tests with small matrices
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=6UL; ++j ) {
            for( size_t k=0UL; k<=i*j; ++k ) {
               for( size_t l=0UL; l<=i*j; ++l ) {
                  RUN_SMATSMATADD_TEST( CMCb( i, j, k ), CMCa( i, j, l ) );
               }
            }
         }
      }

      // Running tests with large matrices
      RUN_SMATSMATADD_TEST( CMCb(  67UL,  67UL,  7UL ), CMCa(  67UL,  67UL,  7UL ) );
      RUN_SMATSMATADD_TEST( CMCb(  67UL, 127UL, 13UL ), CMCa(  67UL, 127UL, 13UL ) );
      RUN_SMATSMATADD_TEST( CMCb( 128UL,  64UL,  8UL ), CMCa( 128UL,  64UL,  8UL ) );
      RUN_SMATSMATADD_TEST( CMCb( 128UL, 128UL, 16UL ), CMCa( 128UL, 128UL, 16UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse matrix/sparse matrix addition:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
