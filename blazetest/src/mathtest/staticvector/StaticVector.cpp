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
   testSubscript();
   testNonZeros();
   testReset();
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
      test_ = "StaticVector default constructor";
   
      blaze::StaticVector<int,5UL,blaze::rowVector> vec;
      
      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 0UL );
   
      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Homogeneous initialization
   {
      test_ = "StaticVector homogeneous initialization constructor";
   
      blaze::StaticVector<int,3UL,blaze::rowVector> vec( 2 );
      
      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );
   
      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // 2D initialization constructor
   {
      test_ = "StaticVector 2D initialization constructor";
   
      blaze::StaticVector<int,2UL,blaze::rowVector> vec( 3, 5 );
      
      checkSize    ( vec, 2UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );
   
      if( vec[0] != 3 || vec[1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // 3D initialization constructor
   {
      test_ = "StaticVector 3D initialization constructor";

      blaze::StaticVector<int,3UL,blaze::rowVector> vec( 3, 5, 2 );
      
      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );
   
      if( vec[0] != 3 || vec[1] != 5 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3 5 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 4D initialization constructor
   {
      test_ = "StaticVector 4D initialization constructor";

      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 3, 5, 2, -7 );
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );
   
      if( vec[0] != 3 || vec[1] != 5 || vec[2] != 2 || vec[3] != -7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3 5 2 -7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 5D initialization constructor
   {
      test_ = "StaticVector 5D initialization constructor";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec( 3, 5, 2, -7, -1 );
      
      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );
   
      if( vec[0] != 3 || vec[1] != 5 || vec[2] != 2 || vec[3] != -7 || vec[4] != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3 5 2 -7 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // 6D initialization constructor
   {
      test_ = "StaticVector 6D initialization constructor";

      blaze::StaticVector<int,6UL,blaze::rowVector> vec( 3, 5, 2, -7, -1, 4 );
      
      checkSize    ( vec, 6UL );
      checkCapacity( vec, 6UL );
      checkNonZeros( vec, 6UL );
   
      if( vec[0] != 3 || vec[1] != 5 || vec[2] != 2 || vec[3] != -7 || vec[4] != -1 || vec[5] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3 5 2 -7 -1 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Array initialization
   {
      test_ = "StaticVector array initialization constructor";

      int array[4] = { 1, 2, 3, 4 };
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( array );
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   // Copy constructor
   {
      test_ = "StaticVector copy constructor";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec1( 1, 2, 3, 4, 5 );
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2( vec1 );
      
      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );
   
      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the StaticVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void StaticVector::testSubscript()
{
   test_ = "StaticVector::operator[]";

   // Writing the first element
   blaze::StaticVector<int,5UL,blaze::rowVector> vec;
   vec[2] = 1;
   
   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 1UL );

   if( vec[2] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Writing the second element
   vec[4] = 2;
   
   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 2UL );

   if( vec[2] != 1 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 2 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Writing the third element
   vec[3] = 3;
   
   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 3UL );
   
   if( vec[2] != 1 || vec[3] != 3 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Writing the fourth element
   vec[0] = 4;
   
   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 1 3 2 )\n";
      throw std::runtime_error( oss.str() );
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
   test_ = "StaticVector::nonZeros()";

   {
      blaze::StaticVector<int,4UL,blaze::rowVector> vec;
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );
   
      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 1, 2, 0, 3 );
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 3UL );
   
      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 3 )\n";
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
   test_ = "StaticVector::reset()";

   // Initialization check
   blaze::StaticVector<int,4UL,blaze::rowVector> vec( 1, 2, 3, 4 );
   
   checkSize    ( vec, 4UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );
   
   if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 3 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   // Resetting the vector
   vec.reset();
   
   checkSize    ( vec, 4UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 0UL );
   
   if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Reset operation failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 0 0 )\n";
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
   test_ = "StaticVector::normalize()";

   // Initialization check
   blaze::StaticVector<double,4UL,blaze::rowVector> vec( 1.0, 2.0, 3.0, 4.0 );
   
   checkSize    ( vec, 4UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );
   
   if( vec[0] != 1.0 || vec[1] != 2.0 || vec[2] != 3.0 || vec[3] != 4.0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 3 4 )\n";
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
      oss << " Test: " << test_ << "\n"
          << " Error: Normalization failed\n"
          << " Details:\n"
          << "   Result: " << vec.length() << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the scale member function of StaticVector.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the scale member function of StaticVector.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void StaticVector::testScale()
{
   test_ = "StaticVector::scale()";

   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec;
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 3;
      vec[3] = 4;
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );
   
      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Integral scaling of the vector
      vec.scale( 2 );
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );
   
      if( vec[0] != 2 || vec[1] != 4 || vec[2] != 6 || vec[3] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 4 6 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Floating point scaling of the vector
      vec.scale( 0.5 );
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );
   
      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
   
   {
      using blaze::complex;
   
      blaze::StaticVector<complex<float>,2UL,blaze::rowVector> vec;
      vec[0] = complex<float>( 1.0F, 0.0F );
      vec[1] = complex<float>( 2.0F, 0.0F );
      vec.scale( complex<float>( 3.0F, 0.0F ) );
      
      checkSize    ( vec, 2UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );
   
      if( vec[0] != complex<float>( 3.0F, 0.0F ) || vec[1] != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( (3,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
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
   test_ = "StaticVector swap";

   blaze::StaticVector<int,4UL,blaze::rowVector> vec1( 1, 2, 3, 4 );
   blaze::StaticVector<int,4UL,blaze::rowVector> vec2( 4, 3, 2, 1 );
   
   swap( vec1, vec2 );
   
   checkSize    ( vec1, 4UL );
   checkCapacity( vec1, 4UL );
   checkNonZeros( vec1, 4UL );
   
   if( vec1[0] != 4 || vec1[1] != 3 || vec1[2] != 2 || vec1[3] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 1 2 3 4 )\n";
      throw std::runtime_error( oss.str() );
   }
   
   checkSize    ( vec2, 4UL );
   checkCapacity( vec2, 4UL );
   checkNonZeros( vec2, 4UL );
   
   if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 4 3 2 1 )\n";
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
   test_ = "min() function";

   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 1, -2, 3, -4 );
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );
   
      if( vec[0] != 1 || vec[1] != -2 || vec[2] != 3 || vec[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 -2 3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the min function
      const int minimum = min( vec );
   
      if( minimum != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
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
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );
   
      if( vec[0] != -1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( -1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the min function
      const int minimum = min( vec );
   
      if( minimum != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
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
   test_ = "max() function";

   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 1, -2, -3, -4 );
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );
   
      if( vec[0] != 1 || vec[1] != -2 || vec[2] != -3 || vec[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 -2 -3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the max function
      const int maximum = max( vec );
   
      if( maximum != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
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
      
      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );
   
      if( vec[0] != -1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( -1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   
      // Testing the max function
      const int maximum = max( vec );
   
      if( maximum != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
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
