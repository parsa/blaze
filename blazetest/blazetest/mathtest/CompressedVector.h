//=================================================================================================
/*!
//  \file blazetest/mathtest/CompressedVector.h
//  \brief Header file for the CompressedVector math test
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

#ifndef _BLAZETEST_MATHTEST_COMPRESSEDVECTOR_H_
#define _BLAZETEST_MATHTEST_COMPRESSEDVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace compressedvector {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the CompressedVector math test.
//
// The CompressedVector class represents a test suite for the blaze::CompressedVector class
// template. It performs a series of both compile time as well as runtime tests.
*/
class CompressedVector
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit CompressedVector();
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
   void testNonZeros    ();
   void testReset       ();
   void testClear       ();
   void testAppend      ();
   void testInsert      ();
   void testErase       ();
   void testFind        ();
   void testResize      ();
   void testReserve     ();
   void testLength      ();
   void testNormalize   ();
   void testScale       ();
   void testSwap        ();
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
   typedef blaze::CompressedVector<int,blaze::rowVector>  VT;   //!< Type of the compressed vector
   typedef typename VT::TransposeType                     TVT;  //!< Transpose compressed vector type
   typedef typename VT::ElementType                       ET;   //!< Element type of the compressed vector
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TVT );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT, typename TVT::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( typename VT::ElementType, typename TVT::ElementType );
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
/*!\brief Checking the size of the given compressed vector.
//
// \param vector The compressed vector to be checked.
// \param expectedSize The expected size of the compressed vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given compressed vector. In case the actual size
// does not correspond to the given expected size, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the compressed vector
void CompressedVector::checkSize( const Type& vector, size_t expectedSize ) const
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
/*!\brief Checking the capacity of the given compressed vector.
//
// \param vector The compressed vector to be checked.
// \param minCapacity The expected minimum capacity of the compressed vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given compressed vector. In case the actual capacity
// is smaller than the given expected minimum capacity, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the compressed vector
void CompressedVector::checkCapacity( const Type& vector, size_t minCapacity ) const
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
/*!\brief Checking the number of non-zero elements of the given compressed vector.
//
// \param vector The compressed vector to be checked.
// \param expectedNonZeros The expected number of non-zero elements of the compressed vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given compressed vector. In
// case the actual number of non-zero elements does not correspond to the given expected number,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the compressed vector
void CompressedVector::checkNonZeros( const Type& vector, size_t expectedNonZeros ) const
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
/*!\brief Testing the functionality of the CompressedVector class template.
//
// \return void
*/
void runTest()
{
   CompressedVector();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the CompressedVector test.
*/
#define RUN_COMPRESSEDVECTOR_TEST \
   blazetest::mathtest::compressedvector::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace compressedvector

} // namespace mathtest

} // namespace blazetest

#endif
