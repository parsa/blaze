//=================================================================================================
/*!
//  \file src/mathtest/inversion/DenseTest.cpp
//  \brief Source file for the dense matrix inversion test
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/inversion/DenseTest.h>


namespace blazetest {

namespace mathtest {

namespace inversion {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest inversion test.
//
// \exception std::runtime_error Matrix inversion error detected.
*/
DenseTest::DenseTest()
{
   using blaze::rowMajor;
   using blaze::columnMajor;

   typedef blaze::complex<float>   cfloat;
   typedef blaze::complex<double>  cdouble;

   testInversion< blaze::StaticMatrix<float,3UL,3UL,rowMajor> >();
   testInversion< blaze::StaticMatrix<float,3UL,3UL,columnMajor> >();
   testInversion< blaze::StaticMatrix<double,3UL,3UL,rowMajor> >();
   testInversion< blaze::StaticMatrix<double,3UL,3UL,columnMajor> >();
   testInversion< blaze::StaticMatrix<cfloat,3UL,3UL,rowMajor> >();
   testInversion< blaze::StaticMatrix<cfloat,3UL,3UL,columnMajor> >();
   testInversion< blaze::StaticMatrix<cdouble,3UL,3UL,rowMajor> >();
   testInversion< blaze::StaticMatrix<cdouble,3UL,3UL,columnMajor> >();

   testInversion< blaze::HybridMatrix<float,3UL,3UL,rowMajor> >();
   testInversion< blaze::HybridMatrix<float,3UL,3UL,columnMajor> >();
   testInversion< blaze::HybridMatrix<double,3UL,3UL,rowMajor> >();
   testInversion< blaze::HybridMatrix<double,3UL,3UL,columnMajor> >();
   testInversion< blaze::HybridMatrix<cfloat,3UL,3UL,rowMajor> >();
   testInversion< blaze::HybridMatrix<cfloat,3UL,3UL,columnMajor> >();
   testInversion< blaze::HybridMatrix<cdouble,3UL,3UL,rowMajor> >();
   testInversion< blaze::HybridMatrix<cdouble,3UL,3UL,columnMajor> >();

   testInversion< blaze::DynamicMatrix<float,rowMajor> >();
   testInversion< blaze::DynamicMatrix<float,columnMajor> >();
   testInversion< blaze::DynamicMatrix<double,rowMajor> >();
   testInversion< blaze::DynamicMatrix<double,columnMajor> >();
   testInversion< blaze::DynamicMatrix<cfloat,rowMajor> >();
   testInversion< blaze::DynamicMatrix<cfloat,columnMajor> >();
   testInversion< blaze::DynamicMatrix<cdouble,rowMajor> >();
   testInversion< blaze::DynamicMatrix<cdouble,columnMajor> >();

   testInversion< blaze::SymmetricMatrix< blaze::DynamicMatrix<float,rowMajor> > >();
   testInversion< blaze::SymmetricMatrix< blaze::DynamicMatrix<float,columnMajor> > >();
   testInversion< blaze::SymmetricMatrix< blaze::DynamicMatrix<double,rowMajor> > >();
   testInversion< blaze::SymmetricMatrix< blaze::DynamicMatrix<double,columnMajor> > >();
   testInversion< blaze::SymmetricMatrix< blaze::DynamicMatrix<cfloat,rowMajor> > >();
   testInversion< blaze::SymmetricMatrix< blaze::DynamicMatrix<cfloat,columnMajor> > >();
   testInversion< blaze::SymmetricMatrix< blaze::DynamicMatrix<cdouble,rowMajor> > >();
   testInversion< blaze::SymmetricMatrix< blaze::DynamicMatrix<cdouble,columnMajor> > >();

   testInversion< blaze::HermitianMatrix< blaze::DynamicMatrix<float,rowMajor> > >();
   testInversion< blaze::HermitianMatrix< blaze::DynamicMatrix<float,columnMajor> > >();
   testInversion< blaze::HermitianMatrix< blaze::DynamicMatrix<double,rowMajor> > >();
   testInversion< blaze::HermitianMatrix< blaze::DynamicMatrix<double,columnMajor> > >();
   testInversion< blaze::HermitianMatrix< blaze::DynamicMatrix<cfloat,rowMajor> > >();
   testInversion< blaze::HermitianMatrix< blaze::DynamicMatrix<cfloat,columnMajor> > >();
   testInversion< blaze::HermitianMatrix< blaze::DynamicMatrix<cdouble,rowMajor> > >();
   testInversion< blaze::HermitianMatrix< blaze::DynamicMatrix<cdouble,columnMajor> > >();

   testInversion< blaze::LowerMatrix< blaze::DynamicMatrix<float,rowMajor> > >();
   testInversion< blaze::LowerMatrix< blaze::DynamicMatrix<float,columnMajor> > >();
   testInversion< blaze::LowerMatrix< blaze::DynamicMatrix<double,rowMajor> > >();
   testInversion< blaze::LowerMatrix< blaze::DynamicMatrix<double,columnMajor> > >();
   testInversion< blaze::LowerMatrix< blaze::DynamicMatrix<cfloat,rowMajor> > >();
   testInversion< blaze::LowerMatrix< blaze::DynamicMatrix<cfloat,columnMajor> > >();
   testInversion< blaze::LowerMatrix< blaze::DynamicMatrix<cdouble,rowMajor> > >();
   testInversion< blaze::LowerMatrix< blaze::DynamicMatrix<cdouble,columnMajor> > >();

   testInversion< blaze::UniLowerMatrix< blaze::DynamicMatrix<float,rowMajor> > >();
   testInversion< blaze::UniLowerMatrix< blaze::DynamicMatrix<float,columnMajor> > >();
   testInversion< blaze::UniLowerMatrix< blaze::DynamicMatrix<double,rowMajor> > >();
   testInversion< blaze::UniLowerMatrix< blaze::DynamicMatrix<double,columnMajor> > >();
   testInversion< blaze::UniLowerMatrix< blaze::DynamicMatrix<cfloat,rowMajor> > >();
   testInversion< blaze::UniLowerMatrix< blaze::DynamicMatrix<cfloat,columnMajor> > >();
   testInversion< blaze::UniLowerMatrix< blaze::DynamicMatrix<cdouble,rowMajor> > >();
   testInversion< blaze::UniLowerMatrix< blaze::DynamicMatrix<cdouble,columnMajor> > >();

   testInversion< blaze::UpperMatrix< blaze::DynamicMatrix<float,rowMajor> > >();
   testInversion< blaze::UpperMatrix< blaze::DynamicMatrix<float,columnMajor> > >();
   testInversion< blaze::UpperMatrix< blaze::DynamicMatrix<double,rowMajor> > >();
   testInversion< blaze::UpperMatrix< blaze::DynamicMatrix<double,columnMajor> > >();
   testInversion< blaze::UpperMatrix< blaze::DynamicMatrix<cfloat,rowMajor> > >();
   testInversion< blaze::UpperMatrix< blaze::DynamicMatrix<cfloat,columnMajor> > >();
   testInversion< blaze::UpperMatrix< blaze::DynamicMatrix<cdouble,rowMajor> > >();
   testInversion< blaze::UpperMatrix< blaze::DynamicMatrix<cdouble,columnMajor> > >();

   testInversion< blaze::UniUpperMatrix< blaze::DynamicMatrix<float,rowMajor> > >();
   testInversion< blaze::UniUpperMatrix< blaze::DynamicMatrix<float,columnMajor> > >();
   testInversion< blaze::UniUpperMatrix< blaze::DynamicMatrix<double,rowMajor> > >();
   testInversion< blaze::UniUpperMatrix< blaze::DynamicMatrix<double,columnMajor> > >();
   testInversion< blaze::UniUpperMatrix< blaze::DynamicMatrix<cfloat,rowMajor> > >();
   testInversion< blaze::UniUpperMatrix< blaze::DynamicMatrix<cfloat,columnMajor> > >();
   testInversion< blaze::UniUpperMatrix< blaze::DynamicMatrix<cdouble,rowMajor> > >();
   testInversion< blaze::UniUpperMatrix< blaze::DynamicMatrix<cdouble,columnMajor> > >();

   testInversion< blaze::DiagonalMatrix< blaze::DynamicMatrix<float,rowMajor> > >();
   testInversion< blaze::DiagonalMatrix< blaze::DynamicMatrix<float,columnMajor> > >();
   testInversion< blaze::DiagonalMatrix< blaze::DynamicMatrix<double,rowMajor> > >();
   testInversion< blaze::DiagonalMatrix< blaze::DynamicMatrix<double,columnMajor> > >();
   testInversion< blaze::DiagonalMatrix< blaze::DynamicMatrix<cfloat,rowMajor> > >();
   testInversion< blaze::DiagonalMatrix< blaze::DynamicMatrix<cfloat,columnMajor> > >();
   testInversion< blaze::DiagonalMatrix< blaze::DynamicMatrix<cdouble,rowMajor> > >();
   testInversion< blaze::DiagonalMatrix< blaze::DynamicMatrix<cdouble,columnMajor> > >();
}
//*************************************************************************************************

} // namespace inversion

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
   std::cout << "   Running dense matrix inversion test..." << std::endl;

   try
   {
      RUN_INVERSION_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix inversion test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
