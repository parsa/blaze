//=================================================================================================
/*!
//  \file src/mathtest/lapack/EigenvalueTest.cpp
//  \brief Source file for the LAPACK eigenvalue test
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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
#include <blazetest/mathtest/lapack/EigenvalueTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace lapack {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the EigenvalueTest class test.
//
// \exception std::runtime_error Eigenvalue computation error detected.
*/
EigenvalueTest::EigenvalueTest()
{
   using blaze::complex;
   using blaze::columnMajor;
   using blaze::rowMajor;


   //=====================================================================================
   // Single precision tests
   //=====================================================================================

   //testGeev < float >();

   testGges < float, columnMajor, columnMajor, columnMajor, columnMajor >();
   testGges < float, columnMajor, columnMajor, columnMajor, rowMajor >();
   testGges < float, columnMajor, columnMajor, rowMajor, columnMajor >();
   testGges < float, columnMajor, columnMajor, rowMajor, rowMajor >();
   testGges < float, columnMajor, rowMajor, columnMajor, columnMajor >();
   testGges < float, columnMajor, rowMajor, columnMajor, rowMajor >();
   testGges < float, columnMajor, rowMajor, rowMajor, columnMajor >();
   testGges < float, columnMajor, rowMajor, rowMajor, rowMajor >();
   testGges < float, rowMajor, columnMajor, columnMajor, columnMajor >();
   testGges < float, rowMajor, columnMajor, columnMajor, rowMajor >();
   testGges < float, rowMajor, columnMajor, rowMajor, columnMajor >();
   testGges < float, rowMajor, columnMajor, rowMajor, rowMajor >();
   testGges < float, rowMajor, rowMajor, columnMajor, columnMajor >();
   testGges < float, rowMajor, rowMajor, columnMajor, rowMajor >();
   testGges < float, rowMajor, rowMajor, rowMajor, columnMajor >();
   testGges < float, rowMajor, rowMajor, rowMajor, rowMajor >();

   testGgesSelect < float, columnMajor, columnMajor, columnMajor, columnMajor >();
   testGgesSelect < float, columnMajor, columnMajor, columnMajor, rowMajor >();
   testGgesSelect < float, columnMajor, columnMajor, rowMajor, columnMajor >();
   testGgesSelect < float, columnMajor, columnMajor, rowMajor, rowMajor >();
   testGgesSelect < float, columnMajor, rowMajor, columnMajor, columnMajor >();
   testGgesSelect < float, columnMajor, rowMajor, columnMajor, rowMajor >();
   testGgesSelect < float, columnMajor, rowMajor, rowMajor, columnMajor >();
   testGgesSelect < float, columnMajor, rowMajor, rowMajor, rowMajor >();
   testGgesSelect < float, rowMajor, columnMajor, columnMajor, columnMajor >();
   testGgesSelect < float, rowMajor, columnMajor, columnMajor, rowMajor >();
   testGgesSelect < float, rowMajor, columnMajor, rowMajor, columnMajor >();
   testGgesSelect < float, rowMajor, columnMajor, rowMajor, rowMajor >();
   testGgesSelect < float, rowMajor, rowMajor, columnMajor, columnMajor >();
   testGgesSelect < float, rowMajor, rowMajor, columnMajor, rowMajor >();
   testGgesSelect < float, rowMajor, rowMajor, rowMajor, columnMajor >();
   testGgesSelect < float, rowMajor, rowMajor, rowMajor, rowMajor >();

   //testSyev < float >();
   //testSyevd< float >();
   //testSyevx< float >();


   //=====================================================================================
   // Double precision tests
   //=====================================================================================

   testGeev < double >();
   
   testGges < double, columnMajor, columnMajor, columnMajor, columnMajor >();
   testGges < double, columnMajor, columnMajor, columnMajor, rowMajor >();
   testGges < double, columnMajor, columnMajor, rowMajor, columnMajor >();
   testGges < double, columnMajor, columnMajor, rowMajor, rowMajor >();
   testGges < double, columnMajor, rowMajor, columnMajor, columnMajor >();
   testGges < double, columnMajor, rowMajor, columnMajor, rowMajor >();
   testGges < double, columnMajor, rowMajor, rowMajor, columnMajor >();
   testGges < double, columnMajor, rowMajor, rowMajor, rowMajor >();
   testGges < double, rowMajor, columnMajor, columnMajor, columnMajor >();
   testGges < double, rowMajor, columnMajor, columnMajor, rowMajor >();
   testGges < double, rowMajor, columnMajor, rowMajor, columnMajor >();
   testGges < double, rowMajor, columnMajor, rowMajor, rowMajor >();
   testGges < double, rowMajor, rowMajor, columnMajor, columnMajor >();
   testGges < double, rowMajor, rowMajor, columnMajor, rowMajor >();
   testGges < double, rowMajor, rowMajor, rowMajor, columnMajor >();
   testGges < double, rowMajor, rowMajor, rowMajor, rowMajor >();

   testGgesSelect < double, columnMajor, columnMajor, columnMajor, columnMajor >();
   testGgesSelect < double, columnMajor, columnMajor, columnMajor, rowMajor >();
   testGgesSelect < double, columnMajor, columnMajor, rowMajor, columnMajor >();
   testGgesSelect < double, columnMajor, columnMajor, rowMajor, rowMajor >();
   testGgesSelect < double, columnMajor, rowMajor, columnMajor, columnMajor >();
   testGgesSelect < double, columnMajor, rowMajor, columnMajor, rowMajor >();
   testGgesSelect < double, columnMajor, rowMajor, rowMajor, columnMajor >();
   testGgesSelect < double, columnMajor, rowMajor, rowMajor, rowMajor >();
   testGgesSelect < double, rowMajor, columnMajor, columnMajor, columnMajor >();
   testGgesSelect < double, rowMajor, columnMajor, columnMajor, rowMajor >();
   testGgesSelect < double, rowMajor, columnMajor, rowMajor, columnMajor >();
   testGgesSelect < double, rowMajor, columnMajor, rowMajor, rowMajor >();
   testGgesSelect < double, rowMajor, rowMajor, columnMajor, columnMajor >();
   testGgesSelect < double, rowMajor, rowMajor, columnMajor, rowMajor >();
   testGgesSelect < double, rowMajor, rowMajor, rowMajor, columnMajor >();
   testGgesSelect < double, rowMajor, rowMajor, rowMajor, rowMajor >();
   
   testSyev < double >();
   testSyevd< double >();
   testSyevx< double >();


   //=====================================================================================
   // Single precision complex tests
   //=====================================================================================

   //testGeev < complex<float> >();
   //testHeev < complex<float> >();
   //testHeevd< complex<float> >();
   //testHeevx< complex<float> >();


   //=====================================================================================
   // Double precision complex tests
   //=====================================================================================

   testGeev < complex<double> >();
   testHeev < complex<double> >();
   testHeevd< complex<double> >();
   testHeevx< complex<double> >();
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
   std::cout << "   Running LAPACK eigenvalue test..." << std::endl;

   try
   {
      RUN_LAPACK_EIGENVALUE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during LAPACK eigenvalue test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
