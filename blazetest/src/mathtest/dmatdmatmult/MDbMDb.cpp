//=================================================================================================
/*!
//  \file src/mathtest/dmatdmatmult/MDbMDb.cpp
//  \brief Source file for the MDbMDb dense matrix/dense matrix multiplication math test
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
#include <blaze/math/DynamicMatrix.h>
#include <blazetest/mathtest/DMatDMatMult.h>
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
   std::cout << "   Running 'MDbMDb'..." << std::endl;

   using blazetest::mathtest::TypeB;

   try
   {
      // Matrix type definitions
      typedef blaze::DynamicMatrix<TypeB>  MDb;

      // Creator type definitions
      typedef blazetest::Creator<MDb>  CMDb;

      // Running tests with small matrices
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=6UL; ++j ) {
            for( size_t k=0UL; k<=6UL; ++k ) {
               RUN_DMATDMATMULT_TEST( CMDb( i, j ), CMDb( j, k ) );
            }
         }
      }

      // Running tests with large matrices
      RUN_DMATDMATMULT_TEST( CMDb(  31UL,  67UL ), CMDb(  67UL,  31UL ) );
      RUN_DMATDMATMULT_TEST( CMDb(  31UL,  67UL ), CMDb(  67UL, 127UL ) );
      RUN_DMATDMATMULT_TEST( CMDb(  67UL,  67UL ), CMDb(  67UL,  67UL ) );
      RUN_DMATDMATMULT_TEST( CMDb( 127UL,  67UL ), CMDb(  67UL,  31UL ) );
      RUN_DMATDMATMULT_TEST( CMDb( 127UL,  67UL ), CMDb(  67UL, 127UL ) );
      RUN_DMATDMATMULT_TEST( CMDb(  32UL,  64UL ), CMDb(  64UL,  32UL ) );
      RUN_DMATDMATMULT_TEST( CMDb(  32UL,  64UL ), CMDb(  64UL, 128UL ) );
      RUN_DMATDMATMULT_TEST( CMDb(  64UL,  64UL ), CMDb(  64UL,  64UL ) );
      RUN_DMATDMATMULT_TEST( CMDb( 128UL,  64UL ), CMDb(  64UL,  32UL ) );
      RUN_DMATDMATMULT_TEST( CMDb( 128UL,  64UL ), CMDb(  64UL, 128UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix/dense matrix multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
