//=================================================================================================
/*!
//  \file blazetest/mathtest/dynamicvector/ClassTest.h
//  \brief Header file for the DynamicVector class test
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

#ifndef _BLAZETEST_MATHTEST_DYNAMICVECTOR_CLASSTEST_H_
#define _BLAZETEST_MATHTEST_DYNAMICVECTOR_CLASSTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/util/AlignmentTrait.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace dynamicvector {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the DynamicVector class template.
//
// This class represents a test suite for the blaze::DynamicVector class template. It performs
// a series of both compile time as well as runtime tests.
*/
class ClassTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit ClassTest();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   template< typename Type >
   void testAlignment( const std::string& type );

   void testConstructors();
   void testAssignment  ();
   void testAddAssign   ();
   void testSubAssign   ();
   void testMultAssign  ();
   void testDivAssign   ();
   void testSubscript   ();
   void testNonZeros    ();
   void testReset       ();
   void testClear       ();
   void testResize      ();
   void testExtend      ();
   void testReserve     ();
   void testScale       ();
   void testSwap        ();
   void testIsDefault   ();
   void testIsNan       ();
   void testLength      ();
   void testNormalize   ();
   void testMinimum     ();
   void testMaximum     ();

   template< typename Type >
   void checkSize( const Type& vector, size_t expectedSize ) const;

   template< typename Type >
   void checkCapacity( const Type& vector, size_t minCapacity ) const;

   template< typename Type >
   void checkNonZeros( const Type& vector, size_t nonzeros ) const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef blaze::DynamicVector<int,blaze::rowVector>  VT;   //!< Type of the dynamic vector
   typedef VT::TransposeType                           TVT;  //!< Transpose dynamic vector type
   typedef VT::ElementType                             ET;   //!< Element type of the dynamic vector
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TVT );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT, TVT::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT::ElementType, TVT::ElementType );
   /*! \endcond */
   //**********************************************************************************************
};
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
// \param type The string representation of the given template type.
// \exception std::runtime_error Error detected.
//
// This function performs a test of the alignment of a DynamicVector instance of the given
// element type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void ClassTest::testAlignment( const std::string& type )
{
   blaze::DynamicVector<Type,blaze::rowVector> vec( 7UL );
   const size_t alignment( blaze::AlignmentTrait<Type>::value );
   const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );

   if( deviation != 0UL ) {
      std::ostringstream oss;
      oss << " Test: DynamicVector" << type << ",rowVector> alignment test\n"
          << " Error: Invalid alignment detected\n"
          << " Details:\n"
          << "   Expected alignment: " << alignment << "\n"
          << "   Deviation         : " << deviation << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the size of the given dynamic vector.
//
// \param vector The dynamic vector to be checked.
// \param expectedSize The expected size of the dynamic vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given dynamic vector. In case the actual size
// does not correspond to the given expected size, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the dynamic vector
void ClassTest::checkSize( const Type& vector, size_t expectedSize ) const
{
   if( vector.size() != expectedSize ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid size detected\n"
          << " Details:\n"
          << "   Size         : " << vector.size() << "\n"
          << "   Expected size: " << expectedSize << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the capacity of the given dynamic vector.
//
// \param vector The dynamic vector to be checked.
// \param minCapacity The expected minimum capacity of the dynamic vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given dynamic vector. In case the actual capacity
// is smaller than the given expected minimum capacity, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the dynamic vector
void ClassTest::checkCapacity( const Type& vector, size_t minCapacity ) const
{
   if( vector.capacity() < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << vector.capacity() << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements of the given dynamic vector.
//
// \param vector The dynamic vector to be checked.
// \param expectedNonZeros The expected number of non-zero elements of the dynamic vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given dynamic vector. In
// case the actual number of non-zero elements does not correspond to the given expected
// number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dynamic vector
void ClassTest::checkNonZeros( const Type& vector, size_t expectedNonZeros ) const
{
   if( vector.nonZeros() != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << vector.nonZeros() << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the functionality of the DynamicVector class template.
//
// \return void
*/
void runTest()
{
   ClassTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the DynamicVector class test.
*/
#define RUN_DYNAMICVECTOR_CLASS_TEST \
   blazetest::mathtest::dynamicvector::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace dynamicvector

} // namespace mathtest

} // namespace blazetest

#endif
