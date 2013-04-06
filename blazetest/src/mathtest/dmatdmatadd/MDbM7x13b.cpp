//=================================================================================================
/*!
//  \file src/mathtest/dmatdmatadd/MDbM7x13b.cpp
//  \brief Source file for the MDbM7x13b dense matrix/dense matrix addition math test
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
#include <blaze/math/StaticMatrix.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/dmatdmatadd/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'MDbM7x13b'..." << std::endl;

   using blazetest::mathtest::TypeB;

   try
   {
      // Matrix type definitions
      typedef blaze::DynamicMatrix<TypeB>          MDb;
      typedef blaze::StaticMatrix<TypeB,7UL,13UL>  M7x13b;

      // Creator type definitions
      typedef blazetest::Creator<MDb>     CMDb;
      typedef blazetest::Creator<M7x13b>  CM7x13b;

      // Running the tests
      RUN_DMATDMATADD_OPERATION_TEST( CMDb( 7UL, 13UL ), CM7x13b() );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix/dense matrix addition:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
