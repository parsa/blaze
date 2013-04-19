//=================================================================================================
/*!
//  \file src/mathtest/smatsvecmult/AliasingTest.cpp
//  \brief Source file for the sparse matrix/sparse vector multiplication aliasing test
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
#include <blazetest/mathtest/smatsvecmult/AliasingTest.h>


namespace blazetest {

namespace mathtest {

namespace smatsvecmult {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the aliasing test class.
//
// \exception std::runtime_error Operation error detected.
*/
AliasingTest::AliasingTest()
   : sA3x4_ ( 3UL, 4UL )
   , sB3x3_ ( 3UL, 3UL )
   , tsA3x4_( 3UL, 4UL )
   , tsB3x3_( 3UL, 3UL )
   , sa4_   ( 4UL )
   , sb4_   ( 4UL )
   , sc3_   ( 3UL )
   , sd3_   ( 3UL )
   , se3_   ( 3UL )
   , da4_   ( 4UL )
   , db3_   ( 3UL )
   , result_()
   , test_  ()
{
   testSMatSVecMult ();
   testTSMatSVecMult();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the sparse matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the sparse matrix/sparse vector multiplication.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testSMatSVecMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "SMatSVecMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = sA3x4_ * sa4_;
      sa4_    = sA3x4_ * sa4_;

      checkResult( sa4_, result_ );
   }

   // Assignment to first operand of left-hand side compound
   {
      test_ = "SMatSVecMult - Assignment to first operand of left-hand side compound";

      initialize();

      result_ = ( sc3_ * trans( sa4_ ) ) * sb4_;
      sc3_    = ( sc3_ * trans( sa4_ ) ) * sb4_;

      checkResult( sc3_, result_ );
   }

   // Assignment to second operand of left-hand side compound
   {
      test_ = "SMatSVecMult - Assignment to second operand of left-hand side compound";

      initialize();

      result_ = ( sc3_ * trans( sa4_ ) ) * sb4_;
      sa4_    = ( sc3_ * trans( sa4_ ) ) * sb4_;

      checkResult( sa4_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "SMatSVecMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = sA3x4_ * ( sa4_ * da4_ );
      sa4_    = sA3x4_ * ( sa4_ * da4_ );

      checkResult( sa4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "SMatSVecMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = sA3x4_ * ( da4_ + sa4_ );
      da4_    = sA3x4_ * ( da4_ + sa4_ );

      checkResult( da4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "SMatSVecMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  sc3_;
      result_ += sB3x3_ * sc3_;
      sc3_    += sB3x3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Addition assignment to first operand of left-hand side compound
   {
      test_ = "SMatSVecMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ += ( sc3_ * trans( sd3_ ) ) * se3_;
      sc3_    += ( sc3_ * trans( sd3_ ) ) * se3_;

      checkResult( sc3_, result_ );
   }

   // Addition assignment to second operand of left-hand side compound
   {
      test_ = "SMatSVecMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sd3_;
      result_ += ( sc3_ * trans( sd3_ ) ) * se3_;
      sd3_    += ( sc3_ * trans( sd3_ ) ) * se3_;

      checkResult( sd3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "SMatSVecMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ += sB3x3_ * ( sc3_ * db3_ );
      sc3_    += sB3x3_ * ( sc3_ * db3_ );

      checkResult( sc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "SMatSVecMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  db3_;
      result_ += sB3x3_ * ( sc3_ * db3_ );
      db3_    += sB3x3_ * ( sc3_ * db3_ );

      checkResult( db3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "SMatSVecMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  sc3_;
      result_ -= sB3x3_ * sc3_;
      sc3_    -= sB3x3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Subtraction assignment to first operand of left-hand side compound
   {
      test_ = "SMatSVecMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ -= ( sc3_ * trans( sd3_ ) ) * se3_;
      sc3_    -= ( sc3_ * trans( sd3_ ) ) * se3_;

      checkResult( sc3_, result_ );
   }

   // Subtraction assignment to second operand of left-hand side compound
   {
      test_ = "SMatSVecMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sd3_;
      result_ -= ( sc3_ * trans( sd3_ ) ) * se3_;
      sd3_    -= ( sc3_ * trans( sd3_ ) ) * se3_;

      checkResult( sd3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "SMatSVecMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ -= sB3x3_ * ( sc3_ * db3_ );
      sc3_    -= sB3x3_ * ( sc3_ * db3_ );

      checkResult( sc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "SMatSVecMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  db3_;
      result_ -= sB3x3_ * ( sc3_ * db3_ );
      db3_    -= sB3x3_ * ( sc3_ * db3_ );

      checkResult( db3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "SMatSVecMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  sc3_;
      result_ *= sB3x3_ * sc3_;
      sc3_    *= sB3x3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Multiplication assignment to first operand of left-hand side compound
   {
      test_ = "SMatSVecMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ *= ( sc3_ * trans( sd3_ ) ) * se3_;
      sc3_    *= ( sc3_ * trans( sd3_ ) ) * se3_;

      checkResult( sc3_, result_ );
   }

   // Multiplication assignment to second operand of left-hand side compound
   {
      test_ = "SMatSVecMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sd3_;
      result_ *= ( sc3_ * trans( sd3_ ) ) * se3_;
      sd3_    *= ( sc3_ * trans( sd3_ ) ) * se3_;

      checkResult( sd3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "SMatSVecMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ *= sB3x3_ * ( sc3_ * db3_ );
      sc3_    *= sB3x3_ * ( sc3_ * db3_ );

      checkResult( sc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "SMatSVecMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  db3_;
      result_ *= sB3x3_ * ( sc3_ * db3_ );
      db3_    *= sB3x3_ * ( sc3_ * db3_ );

      checkResult( db3_, result_ );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the transpose sparse matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the transpose sparse matrix/sparse vector
// multiplication. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
void AliasingTest::testTSMatSVecMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "TSMatSVecMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = tsA3x4_ * sa4_;
      sa4_    = tsA3x4_ * sa4_;

      checkResult( sa4_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "TSMatSVecMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = tsA3x4_ * ( sa4_ * da4_ );
      sa4_    = tsA3x4_ * ( sa4_ * da4_ );

      checkResult( sa4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "TSMatSVecMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = tsA3x4_ * ( da4_ + sa4_ );
      da4_    = tsA3x4_ * ( da4_ + sa4_ );

      checkResult( da4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "TSMatSVecMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  sc3_;
      result_ += tsB3x3_ * sc3_;
      sc3_    += tsB3x3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "TSMatSVecMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ += tsB3x3_ * ( sc3_ * db3_ );
      sc3_    += tsB3x3_ * ( sc3_ * db3_ );

      checkResult( sc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "TSMatSVecMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  db3_;
      result_ += tsB3x3_ * ( sc3_ * db3_ );
      db3_    += tsB3x3_ * ( sc3_ * db3_ );

      checkResult( db3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "TSMatSVecMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  sc3_;
      result_ -= tsB3x3_ * sc3_;
      sc3_    -= tsB3x3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "TSMatSVecMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ -= tsB3x3_ * ( sc3_ * db3_ );
      sc3_    -= tsB3x3_ * ( sc3_ * db3_ );

      checkResult( sc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "TSMatSVecMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  db3_;
      result_ -= tsB3x3_ * ( sc3_ * db3_ );
      db3_    -= tsB3x3_ * ( sc3_ * db3_ );

      checkResult( db3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "TSMatSVecMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  sc3_;
      result_ *= tsB3x3_ * sc3_;
      sc3_    *= tsB3x3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "TSMatSVecMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  sc3_;
      result_ *= tsB3x3_ * ( sc3_ * db3_ );
      sc3_    *= tsB3x3_ * ( sc3_ * db3_ );

      checkResult( sc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "TSMatSVecMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  db3_;
      result_ *= tsB3x3_ * ( sc3_ * db3_ );
      db3_    *= tsB3x3_ * ( sc3_ * db3_ );

      checkResult( db3_, result_ );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member vectors and matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member vectors and matrices to specific predetermined values.
*/
void AliasingTest::initialize()
{
   //=====================================================================================
   // Initialization of the sparse matrices
   //=====================================================================================

   // Initializing the first row-major dense matrix
   sA3x4_(0,0) = -1;
   sA3x4_(0,2) = -2;
   sA3x4_(1,1) =  2;
   sA3x4_(1,2) = -3;
   sA3x4_(1,3) =  1;
   sA3x4_(2,1) =  1;
   sA3x4_(2,2) =  2;
   sA3x4_(2,3) =  2;

   // Initializing the second row-major dense matrix
   sB3x3_(0,0) = -1;
   sB3x3_(1,1) =  1;
   sB3x3_(1,1) = -2;
   sB3x3_(1,1) =  2;
   sB3x3_(2,2) = -3;

   // Initializing the first column-major dense matrix
   tsA3x4_(0,0) = -1;
   tsA3x4_(0,2) = -2;
   tsA3x4_(1,1) =  2;
   tsA3x4_(1,2) = -3;
   tsA3x4_(1,3) =  1;
   tsA3x4_(2,1) =  1;
   tsA3x4_(2,2) =  2;
   tsA3x4_(2,3) =  2;

   // Initializing the second column-major dense matrix
   tsB3x3_(0,0) = -1;
   tsB3x3_(1,1) =  1;
   tsB3x3_(1,1) = -2;
   tsB3x3_(1,1) =  2;
   tsB3x3_(2,2) = -3;


   //=====================================================================================
   // Initialization of the sparse vectors
   //=====================================================================================

   // Initializing the first sparse column vector
   sa4_.resize( 4UL, false );
   sa4_.reset();
   sa4_[0] = -1;
   sa4_[2] = -3;
   sa4_[2] =  2;

   // Initializing the second sparse column vector
   sb4_.resize( 4UL, false );
   sb4_.reset();
   sb4_[1] =  1;
   sb4_[2] =  2;
   sb4_[3] = -1;

   // Initializing the third sparse column vector
   sc3_.resize( 3UL, false );
   sc3_.reset();
   sc3_[0] = 1;
   sc3_[1] = 2;
   sc3_[2] = 3;

   // Initializing the fourth sparse column vector
   sd3_.resize( 3UL, false );
   sd3_.reset();
   sd3_[1] = 2;
   sd3_[2] = 1;

   // Initializing the fifth sparse column vector
   se3_.resize( 3UL, false );
   se3_.reset();
   se3_[1] = 1;
   se3_[2] = 3;


   //=====================================================================================
   // Initialization of the dense vectors
   //=====================================================================================

   // Initializing the first dense column vector
   da4_.resize( 4UL, false );
   da4_[0] = -1;
   da4_[1] =  0;
   da4_[2] = -3;
   da4_[3] =  2;

   // Initializing the second dense column vector
   db3_.resize( 3UL, false );
   db3_[0] = 1;
   db3_[1] = 2;
   db3_[2] = 3;
}
//*************************************************************************************************

} // namespace smatsvecmult

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
   std::cout << "   Running aliasing test..." << std::endl;

   try
   {
      RUN_SMATSVECMULT_ALIASING_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aliasing test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
