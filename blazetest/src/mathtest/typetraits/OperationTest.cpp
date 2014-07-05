//=================================================================================================
/*!
//  \file src/mathtest/typetraits/OperationTest.cpp
//  \brief Source file for the mathematical type traits operation test
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
#include <sstream>
#include <stdexcept>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/typetraits/BaseElementType.h>
#include <blaze/math/typetraits/NumericElementType.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/mathtest/typetraits/OperationTest.h>


namespace blazetest {

namespace mathtest {

namespace typetraits {

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
   testElementType();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST TYPE TRAITS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the mathematical 'BaseElementType' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'BaseElementType' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testBaseElementType()
{
   using blaze::complex;
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::BaseElementType;

   typedef double                                    Type1;  // Built-in data type
   typedef complex<float>                            Type2;  // Complex data type
   typedef StaticVector<int,3UL>                     Type3;  // Vector with built-in element type
   typedef CompressedVector< DynamicVector<float> >  Type4;  // Vector with vector element type

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( BaseElementType<Type1>::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( BaseElementType<Type2>::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( BaseElementType<Type3>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( BaseElementType<Type4>::Type, float );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'NumericElementType' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'NumericElementType' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testNumericElementType()
{
   using blaze::complex;
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::NumericElementType;

   typedef double                                    Type1;  // Built-in data type
   typedef complex<float>                            Type2;  // Complex data type
   typedef StaticVector<int,3UL>                     Type3;  // Vector with built-in element type
   typedef CompressedVector< DynamicVector<float> >  Type4;  // Vector with vector element type

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( NumericElementType<Type1>::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( NumericElementType<Type2>::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( NumericElementType<Type3>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( NumericElementType<Type4>::Type, float );
}
//*************************************************************************************************

} // namespace typetraits

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
   std::cout << "   Running mathematical type traits operation test..." << std::endl;

   try
   {
      RUN_TYPETRAITS_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during mathematical type traits operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
