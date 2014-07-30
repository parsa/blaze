//=================================================================================================
/*!
//  \file src/mathtest/symmetricmatrix/SparseNonNumericTest.cpp
//  \brief Source file for the SymmetricMatrix sparse non-numeric test
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <blaze/math/SparseColumn.h>
#include <blaze/math/SparseRow.h>
#include <blaze/math/SparseSubmatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/symmetricmatrix/SparseNonNumericTest.h>
#include <blazetest/mathtest/IsEqual.h>


namespace blazetest {

namespace mathtest {

namespace symmetricmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SymmetricMatrix sparse non-numeric test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseNonNumericTest::SparseNonNumericTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testScaling();
   testFunctionCall();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testAppend();
   testInsert();
   testErase();
   testResize();
   testReserve();
   testTrim();
   testTranspose();
   testSwap();
   testFind();
   testLowerBound();
   testUpperBound();
   testIsDefault();
   testSubmatrix();
   testRow();
   testColumn();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the SymmetricMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testConstructors()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the SymmetricMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testAssignment()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testAddAssign()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testSubAssign()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testMultAssign()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all SymmetricMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testScaling()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the SymmetricMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseNonNumericTest::testFunctionCall()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testIterator()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testNonZeros()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testReset()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testClear()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testAppend()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testInsert()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testErase()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testResize()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testReserve()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trim() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testTrim()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the SymmetricMatrix
// specialization. Additionally, it performs a test of self-transpose via the \c trans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testTranspose()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the SymmetricMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testSwap()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testFind()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testLowerBound()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testUpperBound()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testIsDefault()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c submatrix() function with the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c submatrix() function with the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testSubmatrix()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c row() function with the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c row() function with the SymmetricMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testRow()
{

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c column() function with the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c column() function with the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseNonNumericTest::testColumn()
{

}
//*************************************************************************************************

} // namespace symmetricmatrix

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
   std::cout << "   Running SymmetricMatrix sparse non-numeric test..." << std::endl;

   try
   {
      RUN_SYMMETRICMATRIX_SPARSENONNUMERIC_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SymmetricMatrix sparse non-numeric test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
