//=================================================================================================
/*!
//  \file src/mathtest/dmatsmatmult/M16x8bMCb.cpp
//  \brief Source file for the M16x8bMCb dense matrix/sparse matrix multiplication math test
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
#include <blazetest/mathtest/DMatSMatMult.h>
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
   std::cout << "   Running 'M16x8bMCb'..." << std::endl;

   using blazetest::mathtest::TypeB;

   try
   {
      // Matrix type definitions
      typedef blaze::StaticMatrix<TypeB,16UL,8UL>  M16x8b;
      typedef blaze::CompressedMatrix<TypeB>       MCb;

      // Creator type definitions
      typedef blazetest::Creator<M16x8b>  CM16x8b;
      typedef blazetest::Creator<MCb>     CMCb;

      // Running the tests
      for( size_t i=0UL; i<=15UL; ++i ) {
         RUN_DMATSMATMULT_TEST( CM16x8b(), CMCb( 8UL, i, 0UL        ) );
         RUN_DMATSMATMULT_TEST( CM16x8b(), CMCb( 8UL, i, 8UL*i*0.25 ) );
         RUN_DMATSMATMULT_TEST( CM16x8b(), CMCb( 8UL, i, 8UL*i*0.5  ) );
         RUN_DMATSMATMULT_TEST( CM16x8b(), CMCb( 8UL, i, 8UL*i*0.75 ) );
         RUN_DMATSMATMULT_TEST( CM16x8b(), CMCb( 8UL, i, 8UL*i      ) );
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
