//=================================================================================================
/*!
//  \file blazetest/mathtest/densecolumn/ClassTest.h
//  \brief Header file for the DenseColumn class test
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

#ifndef _BLAZETEST_MATHTEST_DENSECOLUMN_CLASSTEST_H_
#define _BLAZETEST_MATHTEST_DENSECOLUMN_CLASSTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/views/DenseColumn.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace densecolumn {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the DenseColumn class template.
//
// This class represents a test suite for the blaze::DenseColumn class template. It performs
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
   void testMinimum     ();
   void testMaximum     ();

   template< typename Type >
   void checkSize( const Type& column, size_t expectedSize ) const;

   template< typename Type >
   void checkRows( const Type& matrix, size_t expectedRows ) const;

   template< typename Type >
   void checkColumns( const Type& matrix, size_t expectedColumns ) const;

   template< typename Type >
   void checkCapacity( const Type& object, size_t minCapacity ) const;

   template< typename Type >
   void checkNonZeros( const Type& object, size_t expectedNonZeros ) const;

   template< typename Type >
   void checkNonZeros( const Type& column, size_t index, size_t expectedNonZeros ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initialize();
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  MT;   //!< Row-major dynamic matrix type
   typedef typename MT::OppositeType                  TMT;  //!< Column-major dynamic matrix type
   typedef blaze::DenseColumn<MT>                     CT;   //!< Dense column type for row-major matrices.
   typedef blaze::DenseColumn<TMT>                    TCT;  //!< Dense column type for column-major matrices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT  mat_;   //!< Row-major dynamic matrix.
               /*!< The \f$ 4 \times 5 \f$ matrix is initialized as
                    \f[\left(\begin{array}{*{5}{c}}
                    0 & 0 & -2 &  0 &  7 \\
                    0 & 1 &  0 &  4 & -8 \\
                    0 & 0 & -3 &  5 &  9 \\
                    0 & 0 &  0 & -6 & 10 \\
                    \end{array}\right)\f]. */
   TMT tmat_;  //!< Column-major dynamic matrix.
               /*!< The \f$ 4 \times 5 \f$ matrix is initialized as
                    \f[\left(\begin{array}{*{5}{c}}
                    0 & 0 & -2 &  0 &  7 \\
                    0 & 1 &  0 &  4 & -8 \\
                    0 & 0 & -3 &  5 &  9 \\
                    0 & 0 &  0 & -6 & 10 \\
                    \end{array}\right)\f]. */

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( TMT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( CT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TCT );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( CT  );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( TCT );
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
/*!\brief Checking the size of the given dense column.
//
// \param column The dense column to be checked.
// \param expectedSize The expected size of the dense column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given dense column. In case the actual size does not
// correspond to the given expected size, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense column
void ClassTest::checkSize( const Type& column, size_t expectedSize ) const
{
   if( column.size() != expectedSize ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid size detected\n"
          << " Details:\n"
          << "   Size         : " << column.size() << "\n"
          << "   Expected size: " << expectedSize << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of rows of the given dynamic matrix.
//
// \param matrix The dynamic matrix to be checked.
// \param expectedRows The expected number of rows of the dynamic matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of rows of the given dynamic matrix. In case the actual number
// of rows does not correspond to the given expected number of rows, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >  // Type of the dynamic matrix
void ClassTest::checkRows( const Type& matrix, size_t expectedRows ) const
{
   if( matrix.rows() != expectedRows ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of rows detected\n"
          << " Details:\n"
          << "   Number of rows         : " << matrix.rows() << "\n"
          << "   Expected number of rows: " << expectedRows << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of columns of the given dynamic matrix.
//
// \param matrix The dynamic matrix to be checked.
// \param expectedRows The expected number of columns of the dynamic matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of columns of the given dynamic matrix. In case the
// actual number of columns does not correspond to the given expected number of columns,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dynamic matrix
void ClassTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
{
   if( matrix.columns() != expectedColumns ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of columns detected\n"
          << " Details:\n"
          << "   Number of columns         : " << matrix.columns() << "\n"
          << "   Expected number of columns: " << expectedColumns << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the capacity of the given dense column or dynamic matrix.
//
// \param object The dense column or dynamic matrix to be checked.
// \param minCapacity The expected minimum capacity.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given dense column or dynamic matrix. In
// case the actual capacity is smaller than the given expected minimum capacity, a
// \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense row or dynamic matrix
void ClassTest::checkCapacity( const Type& object, size_t minCapacity ) const
{
   if( object.capacity() < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << object.capacity() << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements of the given dense column or dynamic matrix.
//
// \param object The dense column or dynamic matrix to be checked.
// \param expectedNonZeros The expected number of non-zero elements.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given dense column. In case
// the actual number of non-zero elements does not correspond to the given expected number,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense column or dynamic matrix
void ClassTest::checkNonZeros( const Type& object, size_t expectedNonZeros ) const
{
   if( object.nonZeros() != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << object.nonZeros() << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( object.capacity() < object.nonZeros() ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Number of non-zeros: " << object.nonZeros() << "\n"
          << "   Capacity           : " << object.capacity() << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements in a specific row/column of the given dynamic matrix.
//
// \param matrix The dynamic matrix to be checked.
// \param index The row/column to be checked.
// \param expectedNonZeros The expected number of non-zero elements in the specified row/column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements in the specified row/column of the
// given dynamic matrix. In case the actual number of non-zero elements does not correspond
// to the given expected number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dynamic matrix
void ClassTest::checkNonZeros( const Type& matrix, size_t index, size_t expectedNonZeros ) const
{
   if( matrix.nonZeros( index ) != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements in "
          << ( blaze::IsRowMajorMatrix<Type>::value ? "row " : "column " ) << index << "\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << matrix.nonZeros( index ) << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( matrix.capacity( index ) < matrix.nonZeros( index ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected in "
          << ( blaze::IsRowMajorMatrix<Type>::value ? "row " : "column " ) << index << "\n"
          << " Details:\n"
          << "   Number of non-zeros: " << matrix.nonZeros( index ) << "\n"
          << "   Capacity           : " << matrix.capacity( index ) << "\n";
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
/*!\brief Testing the functionality of the DenseColumn class template.
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
/*!\brief Macro for the execution of the DenseColumn class test.
*/
#define RUN_DENSECOLUMN_CLASS_TEST \
   blazetest::mathtest::densecolumn::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace densecolumn

} // namespace mathtest

} // namespace blazetest

#endif
