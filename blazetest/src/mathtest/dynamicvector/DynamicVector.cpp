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
#include <blaze/math/shims/Equal.h>
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
   testAlignment< signed char          >( "signed char"          );
   testAlignment< unsigned char        >( "unsigned char"        );
   testAlignment< short                >( "short"                );
   testAlignment< unsigned short       >( "unsigned short"       );
   testAlignment< int                  >( "int"                  );
   testAlignment< unsigned int         >( "unsigned int"         );
   testAlignment< float                >( "float"                );
   testAlignment< double               >( "double"               );
   testAlignment< long double          >( "long double"          );
   testAlignment< complex<float>       >( "complex<float>"       );
   testAlignment< complex<double>      >( "complex<double>"      );
   testAlignment< complex<long double> >( "complex<long double>" );

   testConstructors();
   testSubscript();
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
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "DynamicVector default constructor";

      blaze::DynamicVector<int,blaze::rowVector> vec;

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // Size constructor
   //=====================================================================================

   {
      test_ = "DynamicVector size constructor (size 0)";

      blaze::DynamicVector<int,blaze::rowVector> vec( 0UL );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "DynamicVector size constructor (size 10)";

      blaze::DynamicVector<int,blaze::rowVector> vec( 10UL );

      checkSize    ( vec, 10UL );
      checkCapacity( vec, 10UL );
   }


   //=====================================================================================
   // Homogeneous initialization
   //=====================================================================================

   {
      test_ = "DynamicVector homogeneous initialization constructor (size 0)";

      blaze::DynamicVector<int,blaze::rowVector> vec( 0UL, 2 );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "DynamicVector homogeneous initialization constructor (size 3)";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 2 );

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


   //=====================================================================================
   // Array initialization
   //=====================================================================================

   {
      test_ = "DynamicVector array initialization constructor (size 4)";

      int array[4] = { 1, 2, 3, 4 };
      blaze::DynamicVector<int,blaze::rowVector> vec( array );

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


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "DynamicVector copy constructor (size 0)";

      blaze::DynamicVector<int,blaze::rowVector> vec1( 0UL );
      blaze::DynamicVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 0UL );
      checkNonZeros( vec2, 0UL );
   }

   {
      test_ = "DynamicVector copy constructor (size 5)";

      blaze::DynamicVector<int,blaze::rowVector> vec1( 5UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      blaze::DynamicVector<int,blaze::rowVector> vec2( vec1 );

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
/*!\brief Test of the DynamicVector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the DynamicVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DynamicVector::testAssignment()
{
   // Homogeneous assignment
   {
      test_ = "DynamicVector homogeneous assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL );
      vec = 2;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Array assignment
   {
      test_ = "DynamicVector array assignment";

      int array[4] = { 1, 2, 3, 4 };
      blaze::DynamicVector<int,blaze::rowVector> vec;
      vec = array;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Copy assignment
   {
      test_ = "DynamicVector copy assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1( 5UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      blaze::DynamicVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the DynamicVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void DynamicVector::testSubscript()
{
   test_ = "DynamicVector::operator[]";

   // Writing the first element
   blaze::DynamicVector<int,blaze::rowVector> vec( 7UL, 0UL );
   vec[2] = 1;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 1UL );

   if( vec[2] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Writing the second element
   vec[5] = 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 2UL );

   if( vec[2] != 1 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Writing the third element
   vec[3] = 3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 3UL );

   if( vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Writing the fourth element
   vec[0] = 4;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 7UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
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
   test_ = "DynamicVector::nonZeros()";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL, 0 );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] = 1;
      vec[1] = 2;
      vec[2] = 0;
      vec[3] = 3;

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
   test_ = "DynamicVector::reset()";

   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
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
   test_ = "DynamicVector::clear()";

   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
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

   // Clearing the vector
   vec.clear();

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );
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
   test_ = "DynamicVector::resize()";

   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 0
   vec.resize( 0UL );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 3
   vec.resize( 3UL );

   checkSize    ( vec, 3UL );
   checkCapacity( vec, 3UL );

   // Resizing to 5 and preserving the elements
   vec[0] = 1;
   vec[1] = 2;
   vec[2] = 3;
   vec.resize( 5UL, true );

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );

   if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 3 x x )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 2 and preserving the elements
   vec.resize( 2UL, true );

   checkSize    ( vec, 2UL );
   checkCapacity( vec, 2UL );
   checkNonZeros( vec, 2UL );

   if( vec[0] != 1 || vec[1] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 1
   vec.resize( 1UL );

   checkSize    ( vec, 1UL );
   checkCapacity( vec, 1UL );

   // Resizing to 0
   vec.resize( 0 );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );
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
   test_ = "DynamicVector::extend()";

   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Increasing the size of the vector
   vec.extend( 3UL );

   checkSize    ( vec, 3UL );
   checkCapacity( vec, 3UL );

   // Further increasing the size of the vector and preserving the elements
   vec[0] = 1;
   vec[1] = 2;
   vec[2] = 3;
   vec.extend( 2UL, true );

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 3UL );

   if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Extending the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 3 x x )\n";
      throw std::runtime_error( oss.str() );
   }

   // Further increasing the size of the vector
   vec.extend( 10UL, false );

   checkSize    ( vec, 15UL );
   checkCapacity( vec, 15UL );
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
   test_ = "DynamicVector::reserve()";

   // Initialization check
   blaze::DynamicVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Increasing the capacity of the vector
   vec.reserve( 10UL );

   checkSize    ( vec,  0UL );
   checkCapacity( vec, 10UL );
   checkNonZeros( vec,  0UL );

   // Further increasing the capacity of the vector
   vec.reserve( 20UL );

   checkSize    ( vec,  0UL );
   checkCapacity( vec, 20UL );
   checkNonZeros( vec,  0UL );
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
   test_ = "DynamicVector::length()";

   // Initialization check
   blaze::DynamicVector<double,blaze::rowVector> vec( 2UL );
   vec[0] = 3.0;
   vec[1] = 4.0;

   checkSize    ( vec, 2UL );
   checkCapacity( vec, 2UL );
   checkNonZeros( vec, 2UL );

   if( vec[0] != 3.0 || vec[1] != 4.0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 3 4 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Computing the vector length
   const double length( vec.length() );

   if( !blaze::equal( length, 5.0 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
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
   test_ = "DynamicVector::normalize()";

   // Initialization check
   blaze::DynamicVector<double,blaze::rowVector> vec( 4UL );
   vec[0] = 1.0;
   vec[1] = 2.0;
   vec[2] = 3.0;
   vec[3] = 4.0;

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
   test_ = "DynamicVector::scale()";

   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
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

      blaze::DynamicVector<complex<float>,blaze::rowVector> vec( 2UL );
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
   test_ = "DynamicVector swap";

   blaze::DynamicVector<int,blaze::rowVector> vec1( 3UL );
   vec1[0] = 1;
   vec1[1] = 2;
   vec1[2] = 3;

   blaze::DynamicVector<int,blaze::rowVector> vec2( 4UL );
   vec2[0] = 4;
   vec2[1] = 3;
   vec2[2] = 2;
   vec2[3] = 1;

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

   checkSize    ( vec2, 3UL );
   checkCapacity( vec2, 3UL );
   checkNonZeros( vec2, 3UL );

   if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 ) {
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
   test_ = "min() function";

   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  1;
      vec[1] = -2;
      vec[2] =  3;
      vec[3] = -4;

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;

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
   test_ = "max() function";

   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] =  1;
      vec[1] = -2;
      vec[2] = -3;
      vec[3] = -4;

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 4UL );
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;

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
