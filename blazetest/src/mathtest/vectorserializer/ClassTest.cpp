//=================================================================================================
/*!
//  \file src/mathtest/vectorserializer/ClassTest.cpp
//  \brief Source file for the VectorSerializer class test
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

#include <blaze/math/serialization/MatrixSerializer.h>
#include <blaze/math/serialization/VectorSerializer.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/vectorserializer/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace vectorserializer {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the VectorSerializer class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testEmptyVectors();
   testRandomVectors();
   testFailures();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Serialization test with empty vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs serialization tests with empty vectors. In case an error is
// detected, i.e. in case the destination vector is not empty after deserialization, a
// \a std::runtime_error exception is thrown.
*/
void ClassTest::testEmptyVectors()
{
   test_ = "Empty vectors";

   {
      blaze::DynamicVector<int> src;

      runDynamicVectorTests   ( src );
      runCompressedVectorTests( src );
   }

   {
      blaze::DynamicVector<int> src;

      runDynamicVectorTests   ( src );
      runCompressedVectorTests( src );
   }

   {
      blaze::CompressedVector<int> src;

      runDynamicVectorTests   ( src );
      runCompressedVectorTests( src );
   }

   {
      blaze::CompressedVector<int> src;

      runDynamicVectorTests   ( src );
      runCompressedVectorTests( src );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Serialization test with randomly initialized vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs serialization tests with randomly initialized vectors. In
// case an error is detected, i.e. in case a vector cannot be reconstituted from file,
// a \a std::runtime_error exception is thrown.
*/
void ClassTest::testRandomVectors()
{
   test_ = "Randomly initialized vectors";


   //=====================================================================================
   // StaticVector source
   //=====================================================================================

   {
      blaze::StaticVector<int,13UL> src;
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::StaticVector<unsigned int,13UL> src;
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::StaticVector<blaze::complex<float>,13UL> src;
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::StaticVector<blaze::StaticVector<double,3UL>,13UL> src;
      randomize( src );
      runAllTests<13UL>( src );
   }


   //=====================================================================================
   // DynamicVector source
   //=====================================================================================

   {
      blaze::DynamicVector<int> src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::DynamicVector<unsigned int> src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::DynamicVector< blaze::complex<float> > src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::DynamicVector< blaze::StaticVector<double,3UL> > src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }


   //=====================================================================================
   // CompressedVector source
   //=====================================================================================

   {
      blaze::CompressedVector<int> src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::CompressedVector<unsigned int> src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::CompressedVector< blaze::complex<float> > src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::CompressedVector< blaze::StaticVector<double,3UL> > src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of failing serialization attempts.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs tests with failing serialization attempts. In case no error is
// detected, i.e. in case the test is failing, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFailures()
{
   test_ = "Serialization failures";

   try {
      blaze::DynamicMatrix<int> src( 5UL, 1UL );
      blaze::DynamicVector<int> dst;

      randomize( src );

      blaze::Archive<std::stringstream> archive;
      archive << src;
      archive >> dst;

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Type difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicVector<int> src( 5UL );
      blaze::StaticVector<int,3UL> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Size difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicVector<int> src( 5UL );
      blaze::DynamicVector<float> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Element type difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicVector<short> src( 5UL );
      blaze::DynamicVector<long int> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Element size difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}
}
//*************************************************************************************************

} // namespace vectorserializer

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
   std::cout << "   Running VectorSerializer class test..." << std::endl;

   try
   {
      RUN_VECTORSERIALIZER_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during VectorSerializer class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
