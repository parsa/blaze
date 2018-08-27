//=================================================================================================
/*!
//  \file blazetest/utiltest/smallarray/ClassTest.h
//  \brief Header file for the SmallArray test
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZETEST_UTILTEST_SMALLARRAY_CLASSTEST_H_
#define _BLAZETEST_UTILTEST_SMALLARRAY_CLASSTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/util/SmallArray.h>


namespace blazetest {

namespace utiltest {

namespace smallarray {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for the test of the SmallArray class template.
//
// This class represents the collection of tests for the SmallArray class template.
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
   void testConstructors();
   void testAssignment  ();
   void testSubscript   ();
   void testAt          ();
   void testIterator    ();
   void testClear       ();
   void testResize      ();
   void testReserve     ();
   void testShrinkToFit ();
   void testPushBack    ();
   void testInsert      ();
   void testErase       ();
   void testSwap        ();

   template< typename Type >
   void checkSize( const Type& array, size_t expectedSize ) const;

   template< typename Type >
   void checkCapacity( const Type& array, size_t minCapacity ) const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checking the size of the given small array.
//
// \param array The small array to be checked.
// \param expectedSize The expected size of the small array.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given small array. In case the actual size does not
// correspond to the given expected size, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the small array
void ClassTest::checkSize( const Type& array, size_t expectedSize ) const
{
   if( array.size() != expectedSize ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid size detected\n"
          << " Details:\n"
          << "   Size         : " << array.size() << "\n"
          << "   Expected size: " << expectedSize << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the capacity of the given small array.
//
// \param array The small array to be checked.
// \param minCapacity The expected minimum capacity of the small array.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given small array. In case the actual capacity
// is smaller than the given expected minimum capacity, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the small array
void ClassTest::checkCapacity( const Type& array, size_t minCapacity ) const
{
   if( array.capacity() < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << array.capacity() << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for small arrays.
//
// \param os Reference to the output stream.
// \param v Reference to a constant small array object.
// \return Reference to the output stream.
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of elements
inline std::ostream& operator<<( std::ostream& os, const blaze::SmallArray<Type,N>& sv )
{
   os << "(";
   for( const Type& value : sv )
      os << " " << value;
   os << " )";

   return os;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the functionality of the SmallArray class template.
//
// \return void
*/
inline void runTest()
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
/*!\brief Macro for the execution of the SmallArray class test.
*/
#define RUN_SMALLARRAY_CLASS_TEST \
   blazetest::utiltest::smallarray::runTest();
/*! \endcond */
//*************************************************************************************************

} // namespace smallarray

} // namespace utiltest

} // namespace blazetest

#endif
