//=================================================================================================
/*!
//  \file src/mathtest/dmatdvecmult/MDaVDa.cpp
//  \brief Source file for the MDaVDa dense matrix/dense vector multiplication math test
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
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/dmatdvecmult/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'MDaVDa'..." << std::endl;

   using blazetest::mathtest::TypeA;

   try
   {
      // Matrix type definitions
      typedef blaze::DynamicMatrix<TypeA>  MDa;
      typedef blaze::DynamicVector<TypeA>  VDa;

      // Creator type definitions
      typedef blazetest::Creator<MDa>  CMDa;
      typedef blazetest::Creator<VDa>  CVDa;

      // Running tests with small matrices and vectors
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=6UL; ++j ) {
            RUN_DMATDVECMULT_OPERATION_TEST( CMDa( j, i ), CVDa( i ) );
         }
      }

      // Running tests with large matrices and vectors
      RUN_DMATDVECMULT_OPERATION_TEST( CMDa(  67UL, 127UL ), CVDa( 127UL ) );
      RUN_DMATDVECMULT_OPERATION_TEST( CMDa( 127UL,  67UL ), CVDa(  67UL ) );
      RUN_DMATDVECMULT_OPERATION_TEST( CMDa(  64UL, 128UL ), CVDa( 128UL ) );
      RUN_DMATDVECMULT_OPERATION_TEST( CMDa( 128UL,  64UL ), CVDa(  64UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix/dense vector multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
