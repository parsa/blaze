//=================================================================================================
/*!
//  \file src/mathtest/smatdmatadd/MCaM16x8a.cpp
//  \brief Source file for the MCaM16x8a sparse matrix/dense matrix addition math test
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
#include <blazetest/mathtest/smatdmatadd/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'MCaM16x8a'..." << std::endl;

   using blazetest::mathtest::TypeA;

   try
   {
      // Matrix type definitions
      typedef blaze::CompressedMatrix<TypeA>       MCa;
      typedef blaze::StaticMatrix<TypeA,16UL,8UL>  M16x8a;

      // Creator type definitions
      typedef blazetest::Creator<MCa>     CMCa;
      typedef blazetest::Creator<M16x8a>  CM16x8a;

      // Running the tests
      RUN_SMATDMATADD_OPERATION_TEST( CMCa( 16UL, 8UL,   0UL ), CM16x8a() );
      RUN_SMATDMATADD_OPERATION_TEST( CMCa( 16UL, 8UL,  32UL ), CM16x8a() );
      RUN_SMATDMATADD_OPERATION_TEST( CMCa( 16UL, 8UL,  64UL ), CM16x8a() );
      RUN_SMATDMATADD_OPERATION_TEST( CMCa( 16UL, 8UL,  96UL ), CM16x8a() );
      RUN_SMATDMATADD_OPERATION_TEST( CMCa( 16UL, 8UL, 128UL ), CM16x8a() );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse matrix/dense matrix addition:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
