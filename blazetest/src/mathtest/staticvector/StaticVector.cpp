//=================================================================================================
/*!
//  \file src/mathtest/staticvector/StaticVector.cpp
//  \brief Source file for the StaticVector test
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

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/StaticVector.h>
#include <blaze/util/AlignmentTrait.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/StaticVector.h>


namespace blazetest {

namespace mathtest {

namespace staticvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the StaticVector test class.
//
// \exception std::runtime_error Operation error detected.
*/
StaticVector::StaticVector()
{
   testAlignment();
   testConstructors();
   testNonZeros();
   testReset();
   testNormalize();
   testSwap();
   testMinimum();
   testMaximum();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the alignment of different StaticVector instances.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the alignment of different StaticVector instances.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticVector::testAlignment()
{
   // Testing the alignment of a signed integer vector
   {
      blaze::StaticVector<int,7UL,blaze::rowVector> vec;
      const size_t alignment( blaze::AlignmentTrait<int>::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: StaticVector<int> alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of an unsigned integer vector
   {
      blaze::StaticVector<unsigned int,7UL,blaze::rowVector> vec;
      const size_t alignment( blaze::AlignmentTrait<unsigned int>::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: StaticVector<unsigned int> alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of a single precision vector
   {
      blaze::StaticVector<float,7UL,blaze::rowVector> vec;
      const size_t alignment( blaze::AlignmentTrait<float>::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: StaticVector<float> alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of a double precision vector
   {
      blaze::StaticVector<double,7UL,blaze::rowVector> vec;
      const size_t alignment( blaze::AlignmentTrait<double>::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: StaticVector<double> alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of a single precision complex vector
   {
      blaze::StaticVector<blaze::complex<float>,7UL,blaze::rowVector> vec;
      const size_t alignment( blaze::AlignmentTrait< blaze::complex<float> >::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: StaticVector< complex<float> > alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of a double precision complex vector
   {
      blaze::StaticVector<blaze::complex<double>,7UL,blaze::rowVector> vec;
      const size_t alignment( blaze::AlignmentTrait< blaze::complex<double> >::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: StaticVector< complex<double> > alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the StaticVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticVector::testConstructors()
{
   // Default constructor
   {
      blaze::StaticVector<int,5UL,blaze::rowVector> vec;
   
      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector default constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0, 0, 0, 0, 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Homogeneous initialization
   {
      blaze::StaticVector<int,3UL,blaze::rowVector> vec( 2 );
   
      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector homogeneous initialization constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2, 2, 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // 2D initialization constructor
   {
      blaze::StaticVector<int,2UL,blaze::rowVector> vec( 3, 5 );
   
      if( vec[0] != 3 || vec[1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector 2D initialization constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3, 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // 3D initialization constructor
   {
      blaze::StaticVector<int,3UL,blaze::rowVector> vec( 3, 5, 2 );
   
      if( vec[0] != 3 || vec[1] != 5 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector 3D initialization constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3, 5, 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 4D initialization constructor
   {
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 3, 5, 2, -7 );
   
      if( vec[0] != 3 || vec[1] != 5 || vec[2] != 2 || vec[3] != -7 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector 4D initialization constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3, 5, 2, -7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 5D initialization constructor
   {
      blaze::StaticVector<int,5UL,blaze::rowVector> vec( 3, 5, 2, -7, -1 );
   
      if( vec[0] != 3 || vec[1] != 5 || vec[2] != 2 || vec[3] != -7 || vec[4] != -1 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector 5D initialization constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3, 5, 2, -7, -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 6D initialization constructor
   {
      blaze::StaticVector<int,6UL,blaze::rowVector> vec( 3, 5, 2, -7, -1, 4 );
   
      if( vec[0] != 3 || vec[1] != 5 || vec[2] != 2 || vec[3] != -7 || vec[4] != -1 || vec[5] != 4 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector 6D initialization constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3, 5, 2, -7, -1, 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Copy constructor
   {
      blaze::StaticVector<int,5UL,blaze::rowVector> vec1( 1, 2, 3, 4, 5 );
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2( vec1 );
   
      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector copy constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1, 2, 3, 4, 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the nonZeros member function of StaticVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of StaticVector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticVector::testNonZeros()
{
   {
      // Initialization check
      VT vec;
   
      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector::nonZeros()\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0, 0, 0, 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the nonZeros function
      const int nonzeros = vec.nonZeros();
   
      if( nonzeros != 0 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector::nonZeros()\n"
             << " Error: Invalid number of non-zero elements\n"
             << " Details:\n"
             << "   Result: " << nonzeros << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      // Initialization check
      VT vec( 1, 2, 0, 3 );
   
      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector::nonZeros()\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1, 2, 0, 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the nonZeros function
      const int nonzeros = vec.nonZeros();
   
      if( nonzeros != 3 ) {
         std::ostringstream oss;
         oss << " Test: StaticVector::nonZeros()\n"
             << " Error: Invalid number of non-zero elements\n"
             << " Details:\n"
             << "   Result: " << nonzeros << "\n"
             << "   Expected result: 3\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reset member function of StaticVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of StaticVector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticVector::testReset()
{
   // Initialization check
   VT vec( 1, 2, 3, 4 );
   
   if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: StaticVector::reset()\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1, 2, 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Resetting the vector
   vec.reset();
   
   if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: StaticVector::reset()\n"
          << " Error: Reset operation failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0, 0, 0, 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the normalize functionality of the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the normalize and getNormalized functions of the
// StaticVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void StaticVector::testNormalize()
{
   // Initialization check
   blaze::StaticVector<double,4UL,blaze::rowVector> vec( 1.0, 2.0, 3.0, 4.0 );
   
   if( vec[0] != 1.0 || vec[1] != 2.0 || vec[2] != 3.0 || vec[3] != 4.0 ) {
      std::ostringstream oss;
      oss << " Test: StaticVector::getNormalized()\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1, 2, 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Acquiring normalized vector
   const blaze::StaticVector<double,4UL,blaze::rowVector> normalized( vec.getNormalized() );
   
   if( !blaze::equal( normalized.length(), 1.0 ) ) {
      std::ostringstream oss;
      oss << " Test: StaticVector::getNormalized()\n"
          << " Error: Normalization failed\n"
          << " Details:\n"
          << "   Result: " << normalized.length() << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Normalizing the vector
   vec.normalize();
   
   if( !blaze::equal( vec.length(), 1.0 ) ) {
      std::ostringstream oss;
      oss << " Test: StaticVector::normalize()\n"
          << " Error: Normalization failed\n"
          << " Details:\n"
          << "   Result: " << vec.length() << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the swap functionality of the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the swap function of the StaticVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticVector::testSwap()
{
   // Initialization check
   VT vec1( 1, 2, 3, 4 );
   VT vec2( 4, 3, 2, 1 );
   
   if( vec1[0] != 1 || vec1[1] != 2 || vec1[2] != 3 || vec1[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: StaticVector swap\n"
          << " Error: Initialization of first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 1, 2, 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   else if( vec2[0] != 4 || vec2[1] != 3 || vec2[2] != 2 || vec2[3] != 1 ) {
      std::ostringstream oss;
      oss << " Test: StaticVector swap\n"
          << " Error: Initialization of second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 4, 3, 2, 1 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Swapping the vectors
   swap( vec1, vec2 );
   
   if( vec1[0] != 4 || vec1[1] != 3 || vec1[2] != 2 || vec1[3] != 1 ) {
      std::ostringstream oss;
      oss << " Test: StaticVector swap\n"
          << " Error: Swapping the first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 1, 2, 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   else if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: StaticVector swap\n"
          << " Error: Swapping the second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 4, 3, 2, 1 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the min function with the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the min function used with the StaticVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticVector::testMinimum()
{
   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 1, -2, 3, -4 );
   
      if( vec[0] != 1 || vec[1] != -2 || vec[2] != 3 || vec[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: min() function\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1, -2, 3, -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the min function
      const int minimum = min( vec );
   
      if( minimum != -4 ) {
         std::ostringstream oss;
         oss << " Test: min() function\n"
             << " Error: First computation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( -1, 2, 3, 4 );
   
      if( vec[0] != -1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: min() function\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( -1, 2, 3, 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the min function
      const int minimum = min( vec );
   
      if( minimum != -1 ) {
         std::ostringstream oss;
         oss << " Test: min() function\n"
             << " Error: Second computation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: -1\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the max function with the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the max function used with the StaticVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticVector::testMaximum()
{
   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 1, -2, -3, -4 );
   
      if( vec[0] != 1 || vec[1] != -2 || vec[2] != -3 || vec[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: max() function\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1, -2, -3, -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the max function
      const int maximum = max( vec );
   
      if( maximum != 1 ) {
         std::ostringstream oss;
         oss << " Test: max() function\n"
             << " Error: First computation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( -1, 2, 3, 4 );
   
      if( vec[0] != -1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: max() function\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( -1, 2, 3, 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the max function
      const int maximum = max( vec );
   
      if( maximum != 4 ) {
         std::ostringstream oss;
         oss << " Test: max() function\n"
             << " Error: Second computation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace staticvector

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
   std::cout << "   Running StaticVector test..." << std::endl;

   try
   {
      RUN_STATICVECTOR_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during StaticVector test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
