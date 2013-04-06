//=================================================================================================
/*!
//  \file src/mathtest/svectdvecmult/VCaVDb.cpp
//  \brief Source file for the VCaVDb sparse vector/dense vector outer product math test
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
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DynamicVector.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/svectdvecmult/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'VCaVDb'..." << std::endl;

   using blazetest::mathtest::TypeA;
   using blazetest::mathtest::TypeB;

   try
   {
      // Vector type definitions
      typedef blaze::CompressedVector<TypeA>  VCa;
      typedef blaze::DynamicVector<TypeB>     VDb;

      // Creator type definitions
      typedef blazetest::Creator<VCa>  CVCa;
      typedef blazetest::Creator<VDb>  CVDb;

      // Running tests with small vectors
      for( size_t i=0UL; i<=8UL; ++i ) {
         for( size_t j=0UL; j<=8UL; ++j ) {
            for( size_t k=0UL; k<=i; ++k ) {
               RUN_SVECTDVECMULT_OPERATION_TEST( CVCa( i, k ), CVDb( j ) );
            }
         }
      }

      // Running tests with large vectors
      RUN_SVECTDVECMULT_OPERATION_TEST( CVCa(  67UL,  7UL ), CVDb(  67UL ) );
      RUN_SVECTDVECMULT_OPERATION_TEST( CVCa(  67UL, 13UL ), CVDb( 127UL ) );
      RUN_SVECTDVECMULT_OPERATION_TEST( CVCa( 127UL,  7UL ), CVDb(  67UL ) );
      RUN_SVECTDVECMULT_OPERATION_TEST( CVCa( 127UL, 13UL ), CVDb( 127UL ) );
      RUN_SVECTDVECMULT_OPERATION_TEST( CVCa(  64UL,  8UL ), CVDb(  64UL ) );
      RUN_SVECTDVECMULT_OPERATION_TEST( CVCa(  64UL, 16UL ), CVDb( 128UL ) );
      RUN_SVECTDVECMULT_OPERATION_TEST( CVCa( 128UL,  8UL ), CVDb(  64UL ) );
      RUN_SVECTDVECMULT_OPERATION_TEST( CVCa( 128UL, 16UL ), CVDb( 128UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse vector/dense vector outer product:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
