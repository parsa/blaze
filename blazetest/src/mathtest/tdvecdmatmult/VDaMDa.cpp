//=================================================================================================
/*!
//  \file src/mathtest/tdvecdmatmult/VDaMDa.cpp
//  \brief Source file for the VDaMDa dense vector/dense matrix multiplication math test
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
#include <blaze/math/DynamicVector.h>
#include <blazetest/mathtest/TDVecDMatMult.h>
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
   std::cout << "   Running 'VDaMDa'..." << std::endl;

   using blazetest::mathtest::TypeA;

   try
   {
      // Matrix type definitions
      typedef blaze::DynamicVector<TypeA>  VDa;
      typedef blaze::DynamicMatrix<TypeA>  MDa;

      // Creator type definitions
      typedef blazetest::Creator<VDa>  CVDa;
      typedef blazetest::Creator<MDa>  CMDa;

      // Running tests with small vectors and matrices
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=6UL; ++j ) {
            RUN_TDVECDMATMULT_TEST( CVDa( i ), CMDa( i, j ) );
         }
      }

      // Running tests with large vectors and matrices
      RUN_TDVECDMATMULT_TEST( CVDa(  67UL ), CMDa(  67UL, 127UL ) );
      RUN_TDVECDMATMULT_TEST( CVDa( 127UL ), CMDa( 127UL,  67UL ) );
      RUN_TDVECDMATMULT_TEST( CVDa(  64UL ), CMDa(  64UL, 128UL ) );
      RUN_TDVECDMATMULT_TEST( CVDa( 128UL ), CMDa( 128UL,  64UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense vector/dense matrix multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
