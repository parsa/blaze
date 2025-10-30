//=================================================================================================
/*!
//  \file src/mathtest/lapack/DecompositionTest.cpp
//  \brief Source file for the LAPACK decomposition test
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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
#include <blazetest/mathtest/lapack/DecompositionTest.h>

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
/*!\brief Constructor for the DecompositionTest class test.
//
// \exception std::runtime_error Decomposition error detected.
*/
DecompositionTest::DecompositionTest()
{
   using blaze::complex;


   //=====================================================================================
   // Single precision tests
   //=====================================================================================

   //testGetrf< float >();
   //testSytrf< float >();
   //testPotrf< float >();
   //testGeqrf< float >();
   //testOrgqr< float >();
   //testOrg2r< float >();
   //testOrmqr< float >();
   //testGerqf< float >();
   //testOrgrq< float >();
   //testOrgr2< float >();
   //testOrmrq< float >();
   //testGeqlf< float >();
   //testOrgql< float >();
   //testOrg2l< float >();
   //testOrmql< float >();
   //testGelqf< float >();
   //testOrglq< float >();
   //testOrgl2< float >();
   //testOrmlq< float >();


   //=====================================================================================
   // Double precision tests
   //=====================================================================================

   testGetrf< double >();
   testSytrf< double >();
   testPotrf< double >();
   testGeqrf< double >();
   testOrgqr< double >();
   testOrg2r< double >();
   testOrmqr< double >();
   testGerqf< double >();
   testOrgrq< double >();
   testOrgr2< double >();
   testOrmrq< double >();
   testGeqlf< double >();
   testOrgql< double >();
   testOrg2l< double >();
   testOrmql< double >();
   testGelqf< double >();
   testOrglq< double >();
   testOrgl2< double >();
   testOrmlq< double >();


   //=====================================================================================
   // Single precision complex tests
   //=====================================================================================

   //testGetrf< blas_float_complex >();
   //testSytrf< blas_float_complex >();
   //testHetrf< blas_float_complex >();
   //testPotrf< blas_float_complex >();
   //testGeqrf< blas_float_complex >();
   //testUngqr< blas_float_complex >();
   //testUng2r< blas_float_complex >();
   //testUnmqr< blas_float_complex >();
   //testGerqf< blas_float_complex >();
   //testUngrq< blas_float_complex >();
   //testUngr2< blas_float_complex >();
   //testUnmrq< blas_float_complex >();
   //testGeqlf< blas_float_complex >();
   //testUngql< blas_float_complex >();
   //testUng2l< blas_float_complex >();
   //testUnmql< blas_float_complex >();
   //testGelqf< blas_float_complex >();
   //testUnglq< blas_float_complex >();
   //testUngl2< blas_float_complex >();
   //testUnmlq< blas_float_complex >();


   //=====================================================================================
   // Double precision complex tests
   //=====================================================================================

   testGetrf< blas_double_complex >();
   testSytrf< blas_double_complex >();
   testHetrf< blas_double_complex >();
   testPotrf< blas_double_complex >();
   testGeqrf< blas_double_complex >();
   testUngqr< blas_double_complex >();
   testUng2r< blas_double_complex >();
   testUnmqr< blas_double_complex >();
   testGerqf< blas_double_complex >();
   testUngrq< blas_double_complex >();
   testUngr2< blas_double_complex >();
   testUnmrq< blas_double_complex >();
   testGeqlf< blas_double_complex >();
   testUngql< blas_double_complex >();
   testUng2l< blas_double_complex >();
   testUnmql< blas_double_complex >();
   testGelqf< blas_double_complex >();
   testUnglq< blas_double_complex >();
   testUngl2< blas_double_complex >();
   testUnmlq< blas_double_complex >();
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
   std::cout << "   Running LAPACK decomposition test..." << std::endl;

   try
   {
      RUN_LAPACK_DECOMPOSITION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during LAPACK decomposition test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
