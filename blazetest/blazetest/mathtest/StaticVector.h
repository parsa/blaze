//=================================================================================================
/*!
//  \file blazetest/mathtest/StaticVector.h
//  \brief Header file for the StaticVector math test
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

#ifndef _BLAZETEST_MATHTEST_STATICVECTOR_H_
#define _BLAZETEST_MATHTEST_STATICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <typeinfo>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/StaticVector.h>
#include <blaze/util/constraints/SameType.h>


namespace blazetest {

namespace mathtest {

namespace staticvector {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the StaticVector math test.
//
// The StaticVector class represents a test suite for the blaze::StaticVector class template.
// It performs a series of both compile time as well as runtime tests.
*/
class StaticVector
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit StaticVector();
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
   void testNormalize   ();
   void testScale       ();
   void testSwap        ();
   void testMinimum     ();
   void testMaximum     ();
   //@}
   //**********************************************************************************************
   
   //**Type definitions****************************************************************************
   typedef blaze::StaticVector<int,4UL,blaze::rowVector>  VT;   //!< Type of the static vector
   typedef typename VT::TransposeType                     TVT;  //!< Transpose static vector type
   typedef typename VT::ElementType                       ET;   //!< Element type of the static vector
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
/*!\brief Testing the functionality of the StaticVector class template.
//
// \return void
*/
void runTest()
{
   StaticVector();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the StaticVector test.
*/
#define RUN_STATICVECTOR_TEST \
   blazetest::mathtest::staticvector::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace staticvector

} // namespace mathtest

} // namespace blazetest

#endif
