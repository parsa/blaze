//=================================================================================================
/*!
//  \file blazetest/mathtest/densesubvector/AlignedTest.h
//  \brief Header file for the aligned DenseSubvector class test
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

#ifndef _BLAZETEST_MATHTEST_DENSESUBVECTOR_ALIGNEDTEST_H_
#define _BLAZETEST_MATHTEST_DENSESUBVECTOR_ALIGNEDTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/DenseSubvector.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace densesubvector {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the aligned DenseSubvector class template.
//
// This class represents a test suite for the aligned specialization of the blaze::DenseSubvector
// class template. It performs a series of both compile time as well as runtime tests.
*/
class AlignedTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit AlignedTest();
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
   void testConstructors();
   void testAssignment  ();
   void testAddAssign   ();
   void testSubAssign   ();
   void testMultAssign  ();
   void testDivAssign   ();
   void testSubscript   ();
   void testIterator    ();
   void testNonZeros    ();
   void testReset       ();
   void testScale       ();
   void testIsDefault   ();
   void testIsNan       ();
   void testMinimum     ();
   void testMaximum     ();
   void testSubvector   ();

   template< typename Type >
   void checkSize( const Type& vector, size_t expectedSize ) const;

   template< typename Type >
   void checkNonZeros( const Type& vector, size_t expectedNonZeros ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initialize();
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef blaze::DynamicVector<int,blaze::rowVector>  VT;    //!< Dynamic row vector type
   typedef blaze::DenseSubvector<VT,blaze::aligned>    ASVT;  //!< Aligned subvector type for dynamic row vectors.
   typedef blaze::DenseSubvector<VT,blaze::unaligned>  USVT;  //!< Unaligned subvector type for dynamic row vectors.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT vec1_;  //!< First large dynamic row vector.
              /*!< The 64-dimensional dense vector is randomly initialized. */
   VT vec2_;  //!< Second large dynamic row vector.
              /*!< The 64-dimensional dense vector is randomly initialized. */

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ASVT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( USVT );
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
/*!\brief Checking the size of the given dense vector.
//
// \param vector The dense vector to be checked.
// \param expectedSize The expected size of the dense vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given dense vector. In case the actual size does not
// correspond to the given expected size, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense vector
void AlignedTest::checkSize( const Type& vector, size_t expectedSize ) const
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
/*!\brief Checking the number of non-zero elements of the given dense vector.
//
// \param object The dense vector to be checked.
// \param expectedNonZeros The expected number of non-zero elements.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given dense vector. In case
// the actual number of non-zero elements does not correspond to the given expected number,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense vector
void AlignedTest::checkNonZeros( const Type& vector, size_t expectedNonZeros ) const
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

   if( vector.capacity() < vector.nonZeros() ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Number of non-zeros: " << vector.nonZeros() << "\n"
          << "   Capacity           : " << vector.capacity() << "\n";
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
/*!\brief Testing the functionality of the aligned DenseSubvector class template.
//
// \return void
*/
void runTest()
{
   AlignedTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the aligned DenseSubvector class test.
*/
#define RUN_DENSESUBVECTOR_ALIGNED_TEST \
   blazetest::mathtest::densesubvector::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace densesubvector

} // namespace mathtest

} // namespace blazetest

#endif
