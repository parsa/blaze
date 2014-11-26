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
#include <blaze/util/constraints/SameSize.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/typetraits/MakeSigned.h>
#include <blaze/util/typetraits/MakeUnsigned.h>
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
   testMakeSigned();
   testMakeUnsigned();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST TYPE TRAITS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the \c HAS_MEMBER type trait generation macro.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c HAS_MEMBER type trait generation macro.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testHasMember()
{
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

#if !(defined __INTEL_COMPILER) || ( __INTEL_COMPILER >= 1400 )
   BLAZE_STATIC_ASSERT( HasDataType<Type1>::value == 0 );
   BLAZE_STATIC_ASSERT( HasDataType<Type2>::value == 0 );
   BLAZE_STATIC_ASSERT( HasDataType<Type3>::value == 0 );
   BLAZE_STATIC_ASSERT( HasDataType<Type4>::value == 0 );
   BLAZE_STATIC_ASSERT( HasDataType<Type5>::value == 1 );
   BLAZE_STATIC_ASSERT( HasDataType<Type6>::value == 1 );
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c GET_MEMBER type trait generation macro.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c GET_MEMBER type trait generation macro.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testGetMember()
{
#if !(defined __INTEL_COMPILER) || ( __INTEL_COMPILER >= 1400 )
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type5>::Type, Type5::DataType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type6>::Type, Type6::DataType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type7>::Type, int             );
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c MakeSigned type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c MakeSigned type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testMakeSigned()
{
   using blaze::MakeSigned;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<signed char>::Type   , signed char );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<unsigned char>::Type , signed char );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<short>::Type         , short       );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<unsigned short>::Type, short       );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<int>::Type           , int         );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<unsigned int>::Type  , int         );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<long>::Type          , long        );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<unsigned long>::Type , long        );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<const int>::Type         , const int          );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<volatile int>::Type      , volatile int       );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<const volatile int>::Type, const volatile int );

   BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE( MakeSigned<wchar_t>::Type, wchar_t );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c MakeUnsigned type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c MakeUnsigned type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testMakeUnsigned()
{
   using blaze::MakeUnsigned;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<signed char>::Type   , unsigned char  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<unsigned char>::Type , unsigned char  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<short>::Type         , unsigned short );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<unsigned short>::Type, unsigned short );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<int>::Type           , unsigned int   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<unsigned int>::Type  , unsigned int   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<long>::Type          , unsigned long  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<unsigned long>::Type , unsigned long  );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<const int>::Type         , const unsigned int          );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<volatile int>::Type      , volatile unsigned int       );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<const volatile int>::Type, const volatile unsigned int );

   BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE( MakeUnsigned<wchar_t>::Type, wchar_t );
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
