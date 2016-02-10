//=================================================================================================
/*!
//  \file src/mathtest/lapack/OperationTest.cpp
//  \brief Source file for the LAPACK operation test
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
#include <blazetest/mathtest/lapack/OperationTest.h>


namespace blazetest {

namespace mathtest {

namespace lapack {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the OperationTest class test.
//
// \exception std::runtime_error Operation error detected.
*/
OperationTest::OperationTest()
{
   using blaze::complex;


   //=====================================================================================
   // Single precision tests
   //=====================================================================================

   //testGetrf< float >();
   //testSytrf< float >();
   //testPotrf< float >();
   //testGetri< float >();
   //testSytri< float >();
   //testPotri< float >();
   //testTrtri< float >();
   //testGetrs< float >();
   //testSytrs< float >();
   //testPotrs< float >();
   //testTrtrs< float >();
   //testGesv < float >();
   //testSysv < float >();
   //testPosv < float >();
   //testTrsv < float >();
   //testGeqrf< float >();
   //testOrgqr< float >();
   //testGerqf< float >();
   //testOrgrq< float >();


   //=====================================================================================
   // Double precision tests
   //=====================================================================================

   testGetrf< double >();
   testSytrf< double >();
   testPotrf< double >();
   testGetri< double >();
   testSytri< double >();
   testPotri< double >();
   testTrtri< double >();
   testGetrs< double >();
   testSytrs< double >();
   testPotrs< double >();
   testTrtrs< double >();
   testGesv < double >();
   testSysv < double >();
   testPosv < double >();
   testTrsv < double >();
   testGeqrf< double >();
   testOrgqr< double >();
   testGerqf< double >();
   testOrgrq< double >();
   testGelqf< double >();
   testOrglq< double >();


   //=====================================================================================
   // Single precision complex tests
   //=====================================================================================

   //testGetrf< complex<float> >();
   //testSytrf< complex<float> >();
   //testHetrf< complex<float> >();
   //testPotrf< complex<float> >();
   //testGetri< complex<float> >();
   //testSytri< complex<float> >();
   //testHetri< complex<float> >();
   //testPotri< complex<float> >();
   //testTrtri< complex<float> >();
   //testGetrs< complex<float> >();
   //testSytrs< complex<float> >();
   //testHetrs< complex<float> >();
   //testPotrs< complex<float> >();
   //testTrtrs< complex<float> >();
   //testGesv < complex<float> >();
   //testSysv < complex<float> >();
   //testHesv < complex<float> >();
   //testPosv < complex<float> >();
   //testTrsv < complex<float> >();
   //testGeqrf< complex<float> >();
   //testUngqr< complex<float> >();
   //testGerqf< complex<float> >();
   //testUngrq< complex<float> >();


   //=====================================================================================
   // Double precision complex tests
   //=====================================================================================

   testGetrf< complex<double> >();
   testSytrf< complex<double> >();
   testHetrf< complex<double> >();
   testPotrf< complex<double> >();
   testGetri< complex<double> >();
   testSytri< complex<double> >();
   testHetri< complex<double> >();
   testPotri< complex<double> >();
   testTrtri< complex<double> >();
   testGetrs< complex<double> >();
   testSytrs< complex<double> >();
   testHetrs< complex<double> >();
   testPotrs< complex<double> >();
   testTrtrs< complex<double> >();
   testGesv < complex<double> >();
   testSysv < complex<double> >();
   testHesv < complex<double> >();
   testPosv < complex<double> >();
   testTrsv < complex<double> >();
   testGeqrf< complex<double> >();
   testUngqr< complex<double> >();
   testGerqf< complex<double> >();
   testUngrq< complex<double> >();
   testGelqf< complex<double> >();
   testUnglq< complex<double> >();
}
//*************************************************************************************************

} // namespace lapack

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
   std::cout << "   Running LAPACK operation test..." << std::endl;

   try
   {
      RUN_LAPACK_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during LAPACK operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
