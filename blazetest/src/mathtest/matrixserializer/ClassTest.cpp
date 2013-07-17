//=================================================================================================
/*!
//  \file src/mathtest/matrixserializer/ClassTest.cpp
//  \brief Source file for the MatrixSerializer class test
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

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/serialization/MatrixSerializer.h>
#include <blaze/math/serialization/VectorSerializer.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/matrixserializer/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace matrixserializer {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the MatrixSerializer class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testEmptyMatrices();
   testRandomMatrices();
   testFailures();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Serialization test with empty matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs serialization tests with empty matrices. In case an error is
// detected, i.e. in case the destination matrix is not empty after deserialization, a
// \a std::runtime_error exception is thrown.
*/
void ClassTest::testEmptyMatrices()
{
   test_ = "Empty matrices";

   {
      blaze::DynamicMatrix<int> src;

      runDynamicMatrixTests   ( src );
      runCompressedMatrixTests( src );
   }

   {
      blaze::DynamicMatrix<int> src;

      runDynamicMatrixTests   ( src );
      runCompressedMatrixTests( src );
   }

   {
      blaze::CompressedMatrix<int> src;

      runDynamicMatrixTests   ( src );
      runCompressedMatrixTests( src );
   }

   {
      blaze::CompressedMatrix<int> src;

      runDynamicMatrixTests   ( src );
      runCompressedMatrixTests( src );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Serialization test with randomly initialized matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs serialization tests with randomly initialized matrices. In
// case an error is detected, i.e. in case a matrix cannot be reconstituted from file,
// a \a std::runtime_error exception is thrown.
*/
void ClassTest::testRandomMatrices()
{
   test_ = "Randomly initialized matrices";


   //=====================================================================================
   // Row-major StaticMatrix source
   //=====================================================================================

   {
      blaze::StaticMatrix<int,7UL,13UL,blaze::rowMajor> src;
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::StaticMatrix<unsigned int,7UL,13UL,blaze::rowMajor> src;
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::StaticMatrix<blaze::complex<float>,13UL,7UL,blaze::rowMajor> src;
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::StaticMatrix<blaze::StaticVector<double,3UL>,13UL,7UL,blaze::rowMajor> src;
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Column-major StaticMatrix source
   //=====================================================================================

   {
      blaze::StaticMatrix<int,7UL,13UL,blaze::columnMajor> src;
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::StaticMatrix<unsigned int,7UL,13UL,blaze::columnMajor> src;
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::StaticMatrix<blaze::complex<float>,13UL,7UL,blaze::columnMajor> src;
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::StaticMatrix<blaze::StaticVector<double,3UL>,13UL,7UL,blaze::columnMajor> src;
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Row-major DynamicMatrix source
   //=====================================================================================

   {
      blaze::DynamicMatrix<int,blaze::rowMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::DynamicMatrix<unsigned int,blaze::rowMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::DynamicMatrix<blaze::complex<float>,blaze::rowMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::DynamicMatrix<blaze::StaticVector<double,3UL>,blaze::rowMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Column-major DynamicMatrix source
   //=====================================================================================

   {
      blaze::DynamicMatrix<int,blaze::columnMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::DynamicMatrix<unsigned int,blaze::columnMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::DynamicMatrix<blaze::complex<float>,blaze::columnMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::DynamicMatrix<blaze::StaticVector<double,3UL>,blaze::columnMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Row-major CompressedMatrix source
   //=====================================================================================

   {
      blaze::CompressedMatrix<int,blaze::rowMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::CompressedMatrix<unsigned int,blaze::rowMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::CompressedMatrix<blaze::complex<float>,blaze::rowMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::CompressedMatrix<blaze::StaticVector<double,3UL>,blaze::rowMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Column-major CompressedMatrix source
   //=====================================================================================

   {
      blaze::CompressedMatrix<int,blaze::columnMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::CompressedMatrix<unsigned int,blaze::columnMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::CompressedMatrix<blaze::complex<float>,blaze::columnMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::CompressedMatrix<blaze::StaticVector<double,3UL>,blaze::columnMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
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
      blaze::DynamicVector<int> src( 10UL );
      blaze::DynamicMatrix<int> dst;

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
      blaze::DynamicMatrix<int> src( 4UL, 4UL );
      blaze::StaticMatrix<int,3UL,4UL> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Row difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicMatrix<int> src( 3UL, 5UL );
      blaze::StaticMatrix<int,3UL,4UL> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Column difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicMatrix<int> src( 5UL, 4UL );
      blaze::DynamicMatrix<float> dst;

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
      blaze::DynamicMatrix<short> src( 5UL, 4UL );
      blaze::DynamicMatrix<long int> dst;

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

} // namespace matrixserializer

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
   std::cout << "   Running MatrixSerializer class test..." << std::endl;

   try
   {
      RUN_MATRIXSERIALIZER_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during MatrixSerializer class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
