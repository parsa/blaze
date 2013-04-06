//=================================================================================================
/*!
//  \file src/mathtest/dmatdmatmult/MDaMDa.cpp
//  \brief Source file for the MDaMDa dense matrix/dense matrix multiplication math test
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
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/dmatdmatmult/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'MDaMDa'..." << std::endl;

   using blazetest::mathtest::TypeA;

   try
   {
      // Matrix type definitions
      typedef blaze::DynamicMatrix<TypeA>  MDa;

      // Creator type definitions
      typedef blazetest::Creator<MDa>  CMDa;

      // Running tests with small matrices
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=6UL; ++j ) {
            for( size_t k=0UL; k<=6UL; ++k ) {
               RUN_DMATDMATMULT_OPERATION_TEST( CMDa( i, j ), CMDa( j, k ) );
            }
         }
      }

      // Running tests with large matrices
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 15UL, 37UL ), CMDa(  37UL, 15UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 15UL, 37UL ), CMDa(  37UL, 63UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 37UL, 37UL ), CMDa(  37UL, 37UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 63UL, 37UL ), CMDa(  37UL, 15UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 63UL, 37UL ), CMDa(  37UL, 63UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 16UL, 32UL ), CMDa(  32UL, 16UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 16UL, 32UL ), CMDa(  32UL, 64UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 32UL, 32UL ), CMDa(  32UL, 32UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 64UL, 32UL ), CMDa(  32UL, 16UL ) );
      RUN_DMATDMATMULT_OPERATION_TEST( CMDa( 64UL, 32UL ), CMDa(  32UL, 64UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix/dense matrix multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
