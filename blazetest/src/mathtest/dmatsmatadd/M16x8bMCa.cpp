//=================================================================================================
/*!
//  \file src/mathtest/dmatsmatadd/M16x8bMCa.cpp
//  \brief Source file for the M16x8bMCa dense matrix/sparse matrix addition math test
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
#include <blaze/math/StaticMatrix.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/dmatsmatadd/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'M16x8bMCa'..." << std::endl;

   using blazetest::mathtest::TypeA;
   using blazetest::mathtest::TypeB;

   try
   {
      // Matrix type definitions
      typedef blaze::StaticMatrix<TypeB,16UL,8UL>  M16x8b;
      typedef blaze::CompressedMatrix<TypeA>       MCa;

      // Creator type definitions
      typedef blazetest::Creator<M16x8b>  CM16x8b;
      typedef blazetest::Creator<MCa>     CMCa;

      // Running the tests
      RUN_DMATSMATADD_OPERATION_TEST( CM16x8b(), CMCa( 16UL, 8UL,   0UL ) );
      RUN_DMATSMATADD_OPERATION_TEST( CM16x8b(), CMCa( 16UL, 8UL,  32UL ) );
      RUN_DMATSMATADD_OPERATION_TEST( CM16x8b(), CMCa( 16UL, 8UL,  64UL ) );
      RUN_DMATSMATADD_OPERATION_TEST( CM16x8b(), CMCa( 16UL, 8UL,  96UL ) );
      RUN_DMATSMATADD_OPERATION_TEST( CM16x8b(), CMCa( 16UL, 8UL, 128UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix/sparse matrix addition:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
