//=================================================================================================
/*!
//  \file src/mathtest/dmatsmatmult/M3x3bMCa.cpp
//  \brief Source file for the M3x3bMCa dense matrix/sparse matrix multiplication math test
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
#include <blazetest/mathtest/dmatsmatmult/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'M3x3bMCa'..." << std::endl;

   using blazetest::mathtest::TypeA;
   using blazetest::mathtest::TypeB;

   try
   {
      // Matrix type definitions
      typedef blaze::StaticMatrix<TypeB,3UL,3UL>  M3x3b;
      typedef blaze::CompressedMatrix<TypeA>      MCa;

      // Creator type definitions
      typedef blazetest::Creator<M3x3b>  CM3x3b;
      typedef blazetest::Creator<MCa>    CMCa;

      // Running the tests
      for( size_t i=0UL; i<=5UL; ++i ) {
         RUN_DMATSMATMULT_OPERATION_TEST( CM3x3b(), CMCa( 3UL, i, 0UL        ) );
         RUN_DMATSMATMULT_OPERATION_TEST( CM3x3b(), CMCa( 3UL, i, 3UL*i*0.25 ) );
         RUN_DMATSMATMULT_OPERATION_TEST( CM3x3b(), CMCa( 3UL, i, 3UL*i*0.5  ) );
         RUN_DMATSMATMULT_OPERATION_TEST( CM3x3b(), CMCa( 3UL, i, 3UL*i*0.75 ) );
         RUN_DMATSMATMULT_OPERATION_TEST( CM3x3b(), CMCa( 3UL, i, 3UL*i      ) );
      }
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix/sparse matrix multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
