//=================================================================================================
/*!
//  \file src/mathtest/lu/DenseTest.cpp
//  \brief Source file for the dense matrix LU test
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/lu/DenseTest.h>


namespace blazetest {

namespace mathtest {

namespace lu {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest LU test.
//
// \exception std::runtime_error Error during LU decomposition detected.
*/
DenseTest::DenseTest()
{
   typedef blaze::complex<float>   cfloat;
   typedef blaze::complex<double>  cdouble;

   testGeneral<float  >();
   testGeneral<double >();
   testGeneral<cfloat >();
   testGeneral<cdouble>();

   testSymmetric<float  >();
   testSymmetric<double >();
   testSymmetric<cfloat >();
   testSymmetric<cdouble>();

   testHermitian<float  >();
   testHermitian<double >();
   testHermitian<cfloat >();
   testHermitian<cdouble>();

   testLower<float  >();
   testLower<double >();
   testLower<cfloat >();
   testLower<cdouble>();

   testUniLower<float  >();
   testUniLower<double >();
   testUniLower<cfloat >();
   testUniLower<cdouble>();

   testUpper<float  >();
   testUpper<double >();
   testUpper<cfloat >();
   testUpper<cdouble>();

   testUniUpper<float  >();
   testUniUpper<double >();
   testUniUpper<cfloat >();
   testUniUpper<cdouble>();

   testDiagonal<float  >();
   testDiagonal<double >();
   testDiagonal<cfloat >();
   testDiagonal<cdouble>();
}
//*************************************************************************************************

} // namespace lu

} // namespace mathtest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running dense matrix LU decomposition test..." << std::endl;

   try
   {
      RUN_LU_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix LU decomposition test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
