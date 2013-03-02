//=================================================================================================
/*!
//  \file src/mathtest/dynamicvector/DynamicVector.cpp
//  \brief Source file for the DynamicVector test
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
#include <blaze/math/DynamicVector.h>
#include <blaze/util/AlignmentTrait.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/DynamicVector.h>


namespace blazetest {

namespace mathtest {

namespace dynamicvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DynamicVector test class.
//
// \exception std::runtime_error Operation error detected.
*/
DynamicVector::DynamicVector()
{
   testAlignment();
   testConstructors();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testExtend();
   testReserve();
   testLength();
   testNormalize();
   testScale();
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
/*!\brief Test of the alignment of different DynamicVector instances.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the alignment of different DynamicVector instances.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testAlignment()
{
   // Testing the alignment of a signed integer vector
   {
      blaze::DynamicVector<int,blaze::rowVector> vec( 7UL );
      const size_t alignment( blaze::AlignmentTrait<int>::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector<int> alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of an unsigned integer vector
   {
      blaze::DynamicVector<unsigned int,blaze::rowVector> vec( 7UL );
      const size_t alignment( blaze::AlignmentTrait<unsigned int>::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector<unsigned int> alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of a single precision vector
   {
      blaze::DynamicVector<float,blaze::rowVector> vec( 7UL );
      const size_t alignment( blaze::AlignmentTrait<float>::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector<float> alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of a double precision vector
   {
      blaze::DynamicVector<double,blaze::rowVector> vec( 7UL );
      const size_t alignment( blaze::AlignmentTrait<double>::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector<double> alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of a single precision complex vector
   {
      blaze::DynamicVector<blaze::complex<float>,blaze::rowVector> vec( 7UL );
      const size_t alignment( blaze::AlignmentTrait< blaze::complex<float> >::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector< complex<float> > alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation: " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Testing the alignment of a double precision complex vector
   {
      blaze::DynamicVector<blaze::complex<double>,blaze::rowVector> vec( 7UL );
      const size_t alignment( blaze::AlignmentTrait< blaze::complex<double> >::value );
      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );
      
      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector< complex<double> > alignment test\n"
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
/*!\brief Test of the DynamicVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the DynamicVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testConstructors()
{
   // Default constructor
   {
      blaze::DynamicVector<int,blaze::rowVector> vec;
   
      if( vec.size() != 0UL ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector default constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n()\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Homogeneous initialization
   {
      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 2 );
   
      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector homogeneous initialization constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2, 2, 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Array initialization
   {
      int array[4] = { 1, 2, 3, 4 };
      blaze::DynamicVector<int,blaze::rowVector> vec( array );
      
      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector array initialization constructor\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1, 2, 3, 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Copy constructor
   {
      blaze::DynamicVector<int,blaze::rowVector> vec1( 5UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      blaze::DynamicVector<int,blaze::rowVector> vec2( vec1 );
   
      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector copy constructor\n"
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
/*!\brief Test of the nonZeros member function of DynamicVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the nonZeros member function of DynamicVector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testNonZeros()
{
   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );
   
      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector::nonZeros()\n"
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
         oss << " Test: DynamicVector::nonZeros()\n"
             << " Error: Invalid number of non-zero elements\n"
             << " Details:\n"
             << "   Result: " << nonzeros << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 0;
      vec[3] = 3;
   
      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector::nonZeros()\n"
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
         oss << " Test: DynamicVector::nonZeros()\n"
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
/*!\brief Test of the reset member function of DynamicVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reset member function of DynamicVector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testReset()
{
   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
   vec[0] = 1;
   vec[1] = 2;
   vec[2] = 3;
   vec[3] = 4;
   
   if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::reset()\n"
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
      oss << " Test: DynamicVector::reset()\n"
          << " Error: Reset operation failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0, 0, 0, 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the clear member function of DynamicVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the clear member function of DynamicVector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testClear()
{
   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
   vec[0] = 1;
   vec[1] = 2;
   vec[2] = 3;
   vec[3] = 4;
   
   if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::clear()\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1, 2, 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Clearing the vector
   vec.clear();
   
   if( vec.size() != 0UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::clear()\n"
          << " Error: Reset operation failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n()\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the resize member function of DynamicVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the resize member function of DynamicVector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testResize()
{
   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec;
   
   if( vec.size() != 0UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 0\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Increasing the size of the vector
   vec.resize( 3UL );
   
   if( vec.size() != 3UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 3\n";
      throw std::runtime_error( oss.str() );
   }
   if( vec.capacity() < vec.size() ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Invalid vector capacity detected\n"
          << " Details:\n"
          << "   Size    : " << vec.size() << "\n"
          << "   Capacity: " << vec.capacity() << "\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Furth increasing the size of the vector and preserving the elements
   vec[0] = 1;
   vec[1] = 2;
   vec[2] = 3;
   vec.resize( 5UL, true );
   
   if( vec.size() != 5UL || vec[0] != 1 || vec[1] != 2 || vec[2] != 3 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1, 2, 3, x, x )\n";
      throw std::runtime_error( oss.str() );
   }
   if( vec.capacity() < vec.size() ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Invalid vector capacity detected\n"
          << " Details:\n"
          << "   Size    : " << vec.size() << "\n"
          << "   Capacity: " << vec.capacity() << "\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Decreasing the size of the vector and preserving the elements
   vec.resize( 2UL, true );
   
   if( vec.size() != 2UL || vec[0] != 1 || vec[1] != 2 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1, 2 )\n";
      throw std::runtime_error( oss.str() );
   }
   if( vec.capacity() < vec.size() ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Invalid vector capacity detected\n"
          << " Details:\n"
          << "   Size    : " << vec.size() << "\n"
          << "   Capacity: " << vec.capacity() << "\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Further decreasing the size of the vector
   vec.resize( 1UL );
   
   if( vec.size() != 1UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 1\n";
      throw std::runtime_error( oss.str() );
   }
   if( vec.capacity() < vec.size() ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::resize()\n"
          << " Error: Invalid vector capacity detected\n"
          << " Details:\n"
          << "   Size    : " << vec.size() << "\n"
          << "   Capacity: " << vec.capacity() << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the extend member function of DynamicVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the extend member function of DynamicVector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testExtend()
{
   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec;
   
   if( vec.size() != 0UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::extend()\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 0\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Increasing the size of the vector
   vec.extend( 3UL );
   
   if( vec.size() != 3UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::extend()\n"
          << " Error: Extending the vector failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 3\n";
      throw std::runtime_error( oss.str() );
   }
   if( vec.capacity() < vec.size() ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::extend()\n"
          << " Error: Invalid vector capacity detected\n"
          << " Details:\n"
          << "   Size    : " << vec.size() << "\n"
          << "   Capacity: " << vec.capacity() << "\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Furth increasing the size of the vector and preserving the elements
   vec[0] = 1;
   vec[1] = 2;
   vec[2] = 3;
   vec.extend( 2UL, true );
   
   if( vec.size() != 5UL || vec[0] != 1 || vec[1] != 2 || vec[2] != 3 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::extend()\n"
          << " Error: Extending the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1, 2, 3, x, x )\n";
      throw std::runtime_error( oss.str() );
   }
   if( vec.capacity() < vec.size() ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::extend()\n"
          << " Error: Invalid vector capacity detected\n"
          << " Details:\n"
          << "   Size    : " << vec.size() << "\n"
          << "   Capacity: " << vec.capacity() << "\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Furth increasing the size of the vector
   vec.extend( 10UL, true );
   
   if( vec.size() != 15UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::extend()\n"
          << " Error: Extending the vector failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 15\n";
      throw std::runtime_error( oss.str() );
   }
   if( vec.capacity() < vec.size() ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::extend()\n"
          << " Error: Invalid vector capacity detected\n"
          << " Details:\n"
          << "   Size    : " << vec.size() << "\n"
          << "   Capacity: " << vec.capacity() << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the reserve member function of DynamicVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the reserve member function of DynamicVector. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testReserve()
{
   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec;
   
   if( vec.size() != 0UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::reserve()\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 0\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Increasing the capacity of the vector
   vec.reserve( 10UL );
   
   if( vec.size() != 0UL || vec.capacity() < 10UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::reserve()\n"
          << " Error: Reserving capacity for the vector failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 0\n"
          << "   Capacity: " << vec.capacity() << "\n"
          << "   Expected capacity: >= 10\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Further increasing the capacity of the vector
   vec.reserve( 20UL );
   
   if( vec.size() != 0UL || vec.capacity() < 20UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::reserve()\n"
          << " Error: Reserving capacity for the vector failed\n"
          << " Details:\n"
          << "   Size: " << vec.size() << "\n"
          << "   Expected size: 0\n"
          << "   Capacity: " << vec.capacity() << "\n"
          << "   Expected capacity: >= 20\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the length and sqrLength member functions of DynamicVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the length and sqrLength member functions of DynamicVector.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testLength()
{
   // Initialization check
   blaze::DynamicVector<double,blaze::rowVector> vec( 2UL );
   vec[0] = 3.0;
   vec[1] = 4.0;
   
   if( vec[0] != 3.0 || vec[1] != 4.0 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::length()\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Computing the vector length
   const double length( vec.length() );
   
   if( !blaze::equal( length, 5.0 ) ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::length()\n"
          << " Error: Length computation failed\n"
          << " Details:\n"
          << "   Result: " << length << "\n"
          << "   Expected result: 5\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Computing the vector square length
   const double sqrLength( vec.sqrLength() );
   
   if( !blaze::equal( vec.sqrLength(), 25.0 ) ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::sqrLength()\n"
          << " Error: Square length computation failed\n"
          << " Details:\n"
          << "   Result: " << sqrLength << "\n"
          << "   Expected result: 25\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the normalize and getNormalized member functions of the DynamicVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the normalize and getNormalized member functions of
// the DynamicVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void DynamicVector::testNormalize()
{
   // Initialization check
   blaze::DynamicVector<double,blaze::rowVector> vec( 4UL );
   vec[0] = 1.0;
   vec[1] = 2.0;
   vec[2] = 3.0;
   vec[3] = 4.0;
   
   if( vec[0] != 1.0 || vec[1] != 2.0 || vec[2] != 3.0 || vec[3] != 4.0 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::getNormalized()\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1, 2, 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Acquiring normalized vector
   const blaze::DynamicVector<double,blaze::rowVector> normalized( vec.getNormalized() );
   
   if( !blaze::equal( normalized.length(), 1.0 ) ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector::getNormalized()\n"
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
      oss << " Test: DynamicVector::normalize()\n"
          << " Error: Normalization failed\n"
          << " Details:\n"
          << "   Result: " << vec.length() << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of DynamicVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of DynamicVector.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testScale()
{
   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 3;
      vec[3] = 4;
   
      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector::scale()\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1, 2, 3, 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Integral scaling of the vector
      vec.scale( 2 );
   
      if( vec[0] != 2 || vec[1] != 4 || vec[2] != 6 || vec[3] != 8 ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector::scale()\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2, 4, 6, 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Floating point scaling of the vector
      vec.scale( 0.5 );
   
      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector::scale()\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1, 2, 3, 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   {
      using blaze::complex;
   
      blaze::DynamicVector<complex<float>,blaze::rowVector> vec( 2UL );
      vec[0] = complex<float>( 1.0F, 0.0F );
      vec[1] = complex<float>( 2.0F, 0.0F );
      vec.scale( complex<float>( 3.0F, 0.0F ) );
   
      if( vec[0] != complex<float>( 3.0F, 0.0F ) || vec[1] != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: DynamicVector::scale()\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( (3,0), (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the swap functionality of the DynamicVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the swap function of the DynamicVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testSwap()
{
   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec1( 4UL );
   vec1[0] = 1;
   vec1[1] = 2;
   vec1[2] = 3;
   vec1[3] = 4;

   blaze::DynamicVector<int,blaze::rowVector> vec2( 4UL );
   vec2[0] = 4;
   vec2[1] = 3;
   vec2[2] = 2;
   vec2[3] = 1;
   
   if( vec1[0] != 1 || vec1[1] != 2 || vec1[2] != 3 || vec1[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector swap\n"
          << " Error: Initialization of first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 1, 2, 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   else if( vec2[0] != 4 || vec2[1] != 3 || vec2[2] != 2 || vec2[3] != 1 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector swap\n"
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
      oss << " Test: DynamicVector swap\n"
          << " Error: Swapping the first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 1, 2, 3, 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   else if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector swap\n"
          << " Error: Swapping the second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 4, 3, 2, 1 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the min function with the DynamicVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the min function used with the DynamicVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testMinimum()
{
   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  1;
      vec[1] = -2;
      vec[2] =  3;
      vec[3] = -4;
   
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
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;
   
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
/*!\brief Test of the max function with the DynamicVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the max function used with the DynamicVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testMaximum()
{
   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  1;
      vec[1] = -2;
      vec[2] = -3;
      vec[3] = -4;
   
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
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;
   
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

} // namespace dynamicvector

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
   std::cout << "   Running DynamicVector test..." << std::endl;

   try
   {
      RUN_DYNAMICVECTOR_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during DynamicVector test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
