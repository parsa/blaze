//=================================================================================================
/*!
//  \file src/mathtest/dvecdvecadd/V2bV2a.cpp
//  \brief Source file for the V2bV2a dense vector/dense vector addition math test
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
#include <blaze/math/StaticVector.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/dvecdvecadd/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'V2bV2a'..." << std::endl;

   using blazetest::mathtest::TypeA;
   using blazetest::mathtest::TypeB;

   try
   {
      // Vector type definitions
      typedef blaze::StaticVector<TypeB,2UL>  V2b;
      typedef blaze::StaticVector<TypeA,2UL>  V2a;

      // Creator type definitions
      typedef blazetest::Creator<V2b>  CV2b;
      typedef blazetest::Creator<V2a>  CV2a;

      // Running the tests
      RUN_DVECDVECADD_OPERATION_TEST( CV2b(), CV2a() );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense vector/dense vector addition:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
