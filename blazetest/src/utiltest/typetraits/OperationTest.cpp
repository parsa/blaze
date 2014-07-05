//=================================================================================================
/*!
//  \file src/utiltest/typetraits/OperationTest.cpp
//  \brief Source file for the type traits operation test
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
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/StaticAssert.h>
#include <blazetest/utiltest/typetraits/OperationTest.h>


namespace blazetest {

namespace utiltest {

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
   testHasMember();
   testGetMember();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST TYPE TRAITS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'HAS_MEMBER' type trait generation macro.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'HAS_MEMBER' type trait generation macro. In
// case an error is detected, a compilation error is created.
*/
void OperationTest::testHasMember()
{
   struct Type1 { int value_; };
   class  Type2 { volatile int value_; };
   struct Type3 { void compute(); };
   class  Type4 { void compute() const; };
   struct Type5 { typedef int  DataType; };
   class  Type6 { typedef const int  DataType; };

   BLAZE_STATIC_ASSERT( HasValue<Type1>::value == 1 );
   BLAZE_STATIC_ASSERT( HasValue<Type2>::value == 1 );
   BLAZE_STATIC_ASSERT( HasValue<Type3>::value == 0 );
   BLAZE_STATIC_ASSERT( HasValue<Type4>::value == 0 );
   BLAZE_STATIC_ASSERT( HasValue<Type5>::value == 0 );
   BLAZE_STATIC_ASSERT( HasValue<Type6>::value == 0 );

   BLAZE_STATIC_ASSERT( HasCompute<Type1>::value == 0 );
   BLAZE_STATIC_ASSERT( HasCompute<Type2>::value == 0 );
   BLAZE_STATIC_ASSERT( HasCompute<Type3>::value == 1 );
   BLAZE_STATIC_ASSERT( HasCompute<Type4>::value == 1 );
   BLAZE_STATIC_ASSERT( HasCompute<Type5>::value == 0 );
   BLAZE_STATIC_ASSERT( HasCompute<Type6>::value == 0 );

   BLAZE_STATIC_ASSERT( HasDataType<Type1>::value == 0 );
   BLAZE_STATIC_ASSERT( HasDataType<Type2>::value == 0 );
   BLAZE_STATIC_ASSERT( HasDataType<Type3>::value == 0 );
   BLAZE_STATIC_ASSERT( HasDataType<Type4>::value == 0 );
   BLAZE_STATIC_ASSERT( HasDataType<Type5>::value == 1 );
   BLAZE_STATIC_ASSERT( HasDataType<Type6>::value == 1 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'GET_MEMBER' type trait generation macro.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'GET_MEMBER' type trait generation macro. In
// case an error is detected, a compilation error is created.
*/
void OperationTest::testGetMember()
{
   struct Type1 { typedef float  DataType; };
   struct Type2 { typedef const double  DataType; };
   struct Type3 {};

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type1>::Type, Type1::DataType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type2>::Type, Type2::DataType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type3>::Type, int             );
}
//*************************************************************************************************

} // namespace typetraits

} // namespace utiltest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running type traits operation test..." << std::endl;

   try
   {
      RUN_TYPETRAITS_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during type traits operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
