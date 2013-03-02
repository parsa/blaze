//=================================================================================================
/*!
//  \file blazetest/mathtest/DynamicVector.h
//  \brief Header file for the DynamicVector math test
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

#ifndef _BLAZETEST_MATHTEST_DYNAMICVECTOR_H_
#define _BLAZETEST_MATHTEST_DYNAMICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <typeinfo>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/util/constraints/SameType.h>


namespace blazetest {

namespace mathtest {

namespace dynamicvector {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the DynamicVector math test.
//
// The DynamicVector class represents a test suite for the blaze::DynamicVector class template.
// It performs a series of both compile time as well as runtime tests.
*/
class DynamicVector
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit DynamicVector();
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
   void testAlignment   ();
   void testConstructors();
   void testNonZeros    ();
   void testReset       ();
   void testClear       ();
   void testResize      ();
   void testExtend      ();
   void testReserve     ();
   void testLength      ();
   void testNormalize   ();
   void testScale       ();
   void testSwap        ();
   void testMinimum     ();
   void testMaximum     ();
   //@}
   //**********************************************************************************************
   
   //**Type definitions****************************************************************************
   typedef blaze::DynamicVector<int,blaze::rowVector>  VT;   //!< Type of the dynamic vector
   typedef typename VT::TransposeType                  TVT;  //!< Transpose dynamic vector type
   typedef typename VT::ElementType                    ET;   //!< Element type of the dynamic vector
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TVT );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT, typename TVT::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( typename VT::ElementType, typename TVT::ElementType );
   /*! \endcond */
   //**********************************************************************************************
};
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
   DynamicVector();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the DynamicVector test.
*/
#define RUN_DYNAMICVECTOR_TEST \
   blazetest::mathtest::dynamicvector::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace dynamicvector

} // namespace mathtest

} // namespace blazetest

#endif
