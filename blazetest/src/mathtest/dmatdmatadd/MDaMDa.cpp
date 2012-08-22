//=================================================================================================
/*!
//  \file src/mathtest/dmatdmatadd/MDaMDa.cpp
//  \brief Source file for the MDaMDa dense matrix/dense matrix addition math test
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
#include <blazetest/mathtest/DMatDMatAdd.h>
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
   std::cout << "   Running 'MDaMDa'..." << std::endl;

   using blazetest::mathtest::TypeA;

   try
   {
      // Matrix type definitions
      typedef blaze::DynamicMatrix<TypeA>  MDa;

      // Creator type definitions
      typedef blazetest::Creator<MDa>  CMDa;

      // Running tests with small matrices
      for( size_t i=0UL; i<=9UL; ++i ) {
         for( size_t j=0UL; j<=9UL; ++j ) {
            RUN_DMATDMATADD_TEST( CMDa( i, j ), CMDa( i, j ) );
         }
      }

      // Running tests with large matrices
      RUN_DMATDMATADD_TEST( CMDa(  67UL,  67UL ), CMDa(  67UL,  67UL ) );
      RUN_DMATDMATADD_TEST( CMDa(  67UL, 127UL ), CMDa(  67UL, 127UL ) );
      RUN_DMATDMATADD_TEST( CMDa( 128UL,  64UL ), CMDa( 128UL,  64UL ) );
      RUN_DMATDMATADD_TEST( CMDa( 128UL, 128UL ), CMDa( 128UL, 128UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix/dense matrix addition:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
