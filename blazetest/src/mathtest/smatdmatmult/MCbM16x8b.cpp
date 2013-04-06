//=================================================================================================
/*!
//  \file src/mathtest/smatdmatmult/MCbM16x8b.cpp
//  \brief Source file for the MCbM16x8b sparse matrix/dense matrix multiplication math test
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
#include <blazetest/mathtest/smatdmatmult/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'MCbM16x8b'..." << std::endl;

   using blazetest::mathtest::TypeB;

   try
   {
      // Matrix type definitions
      typedef blaze::CompressedMatrix<TypeB>       MCb;
      typedef blaze::StaticMatrix<TypeB,16UL,8UL>  M16x8b;

      // Creator type definitions
      typedef blazetest::Creator<MCb>     CMCb;
      typedef blazetest::Creator<M16x8b>  CM16x8b;

      // Running the tests
      for( size_t i=0UL; i<=12UL; ++i ) {
         RUN_SMATDMATMULT_OPERATION_TEST( CMCb( i, 16UL, 0UL         ), CM16x8b() );
         RUN_SMATDMATMULT_OPERATION_TEST( CMCb( i, 16UL, i*16UL*0.25 ), CM16x8b() );
         RUN_SMATDMATMULT_OPERATION_TEST( CMCb( i, 16UL, i*16UL*0.5  ), CM16x8b() );
         RUN_SMATDMATMULT_OPERATION_TEST( CMCb( i, 16UL, i*16UL*0.75 ), CM16x8b() );
         RUN_SMATDMATMULT_OPERATION_TEST( CMCb( i, 16UL, i*16UL      ), CM16x8b() );
      }
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse matrix/dense matrix multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
