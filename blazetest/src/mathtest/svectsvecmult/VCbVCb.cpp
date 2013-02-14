//=================================================================================================
/*!
//  \file src/mathtest/svectsvecmult/VCbVCb.cpp
//  \brief Source file for the VCbVCb sparse vector/sparse vector outer product math test
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
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/SVecTSVecMult.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'VCbVCb'..." << std::endl;

   using blazetest::mathtest::TypeB;

   try
   {
      // Vector type definitions
      typedef blaze::CompressedVector<TypeB>  VCb;

      // Creator type definitions
      typedef blazetest::Creator<VCb>  CVCb;

      // Running tests with small vectors
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=6UL; ++j ) {
            for( size_t k=0UL; k<=i; ++k ) {
               for( size_t l=0UL; l<=j; ++l ) {
                  RUN_SVECTSVECMULT_TEST( CVCb( i, k ), CVCb( j, l ) );
               }
            }
         }
      }

      // Running tests with large vectors
      RUN_SVECTSVECMULT_TEST( CVCb(  67UL,  7UL ), CVCb(  67UL,  7UL ) );
      RUN_SVECTSVECMULT_TEST( CVCb(  67UL, 13UL ), CVCb( 127UL, 13UL ) );
      RUN_SVECTSVECMULT_TEST( CVCb( 127UL,  7UL ), CVCb(  67UL,  7UL ) );
      RUN_SVECTSVECMULT_TEST( CVCb( 127UL, 13UL ), CVCb( 127UL, 13UL ) );
      RUN_SVECTSVECMULT_TEST( CVCb(  64UL,  8UL ), CVCb(  64UL,  8UL ) );
      RUN_SVECTSVECMULT_TEST( CVCb(  64UL, 16UL ), CVCb( 128UL, 16UL ) );
      RUN_SVECTSVECMULT_TEST( CVCb( 128UL,  8UL ), CVCb(  64UL,  8UL ) );
      RUN_SVECTSVECMULT_TEST( CVCb( 128UL, 16UL ), CVCb( 128UL, 16UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse vector/sparse vector outer product:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
