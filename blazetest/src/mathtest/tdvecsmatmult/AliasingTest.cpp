//=================================================================================================
/*!
//  \file src/mathtest/tdvecsmatmult/AliasingTest.cpp
//  \brief Source file for the dense vector/sparse matrix multiplication aliasing test
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
#include <blazetest/mathtest/tdvecsmatmult/AliasingTest.h>


namespace blazetest {

namespace mathtest {

namespace tdvecsmatmult {

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
   : sA4x3_ ( 4UL, 3UL )
   , sB3x3_ ( 3UL, 3UL )
   , tsA4x3_( 4UL, 3UL )
   , tsB3x3_( 3UL, 3UL )
   , tda4_  ( 4UL )
   , tdb4_  ( 4UL )
   , tdc3_  ( 3UL )
   , tdd3_  ( 3UL )
   , tsa4_  ( 4UL )
   , tsb3_  ( 3UL )
   , result_()
   , test_  ()
{
   testTDVecSMatMult ();
   testTDVecTSMatMult();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the dense vector/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the dense vector/sparse matrix multiplication.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testTDVecSMatMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "TDVecSMatMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = tda4_ * sA4x3_;
      tda4_   = tda4_ * sA4x3_;

      checkResult( tda4_, result_ );
   }

   // Assignment to first operand of left-hand side compound
   {
      test_ = "TDVecSMatMult - Assignment to first operand of left-hand side compound";

      initialize();

      result_ = tdb4_ * ( trans( tda4_ ) * tsb3_ );
      tda4_   = tdb4_ * ( trans( tda4_ ) * tsb3_ );

      checkResult( tda4_, result_ );
   }

   // Assignment to second operand of left-hand side compound
   {
      test_ = "TDVecSMatMult - Assignment to second operand of left-hand side compound";

      initialize();

      result_ = tdb4_ * ( trans( tda4_ ) * tsb3_ );
      tsb3_   = tdb4_ * ( trans( tda4_ ) * tsb3_ );

      checkResult( tsb3_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "TDVecSMatMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = ( tda4_ + tsa4_ ) * sA4x3_;
      tda4_   = ( tda4_ + tsa4_ ) * sA4x3_;

      checkResult( tda4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "TDVecSMatMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = ( tda4_ + tsa4_ ) * sA4x3_;
      tsa4_   = ( tda4_ + tsa4_ ) * sA4x3_;

      checkResult( tsa4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "TDVecSMatMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ += tdc3_ * sB3x3_;
      tdc3_   += tdc3_ * sB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to first operand of left-hand side compound
   {
      test_ = "TDVecSMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ += tdd3_ * ( trans( tdc3_ ) * tsb3_ );
      tdc3_   += tdd3_ * ( trans( tdc3_ ) * tsb3_ );

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to second operand of left-hand side compound
   {
      test_ = "TDVecSMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ += tdd3_ * ( trans( tdc3_ ) * tsb3_ );
      tsb3_   += tdd3_ * ( trans( tdc3_ ) * tsb3_ );

      checkResult( tsb3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "TDVecSMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ += ( tdc3_ + tsb3_ ) * sB3x3_;
      tdc3_   += ( tdc3_ + tsb3_ ) * sB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "TDVecSMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ += ( tdc3_ + tsb3_ ) * sB3x3_;
      tsb3_   += ( tdc3_ + tsb3_ ) * sB3x3_;

      checkResult( tsb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "TDVecSMatMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ -= tdc3_ * sB3x3_;
      tdc3_   -= tdc3_ * sB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to first operand of left-hand side compound
   {
      test_ = "TDVecSMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ -= tdd3_ * ( trans( tdc3_ ) * tsb3_ );
      tdc3_   -= tdd3_ * ( trans( tdc3_ ) * tsb3_ );

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to second operand of left-hand side compound
   {
      test_ = "TDVecSMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ -= tdd3_ * ( trans( tdc3_ ) * tsb3_ );
      tsb3_   -= tdd3_ * ( trans( tdc3_ ) * tsb3_ );

      checkResult( tsb3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "TDVecSMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ -= ( tdc3_ + tsb3_ ) * sB3x3_;
      tdc3_   -= ( tdc3_ + tsb3_ ) * sB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "TDVecSMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ -= ( tdc3_ + tsb3_ ) * sB3x3_;
      tsb3_   -= ( tdc3_ + tsb3_ ) * sB3x3_;

      checkResult( tsb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "TDVecSMatMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ *= tdc3_ * sB3x3_;
      tdc3_   *= tdc3_ * sB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to first operand of left-hand side compound
   {
      test_ = "TDVecSMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ *= tdd3_ * ( trans( tdc3_ ) * tsb3_ );
      tdc3_   *= tdd3_ * ( trans( tdc3_ ) * tsb3_ );

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to second operand of left-hand side compound
   {
      test_ = "TDVecSMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ *= tdd3_ * ( trans( tdc3_ ) * tsb3_ );
      tsb3_   *= tdd3_ * ( trans( tdc3_ ) * tsb3_ );

      checkResult( tsb3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "TDVecSMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ *= ( tdc3_ + tsb3_ ) * sB3x3_;
      tdc3_   *= ( tdc3_ + tsb3_ ) * sB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "TDVecSMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ *= ( tdc3_ + tsb3_ ) * sB3x3_;
      tsb3_   *= ( tdc3_ + tsb3_ ) * sB3x3_;

      checkResult( tsb3_, result_ );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the transpose dense vector/transpose dense matrix multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the transpose dense vector/transpose dense matrix
// multiplication. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testTDVecTSMatMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "TDVecTSMatMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = tda4_ * tsA4x3_;
      tda4_   = tda4_ * tsA4x3_;

      checkResult( tda4_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "TDVecTSMatMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = ( tda4_ + tsa4_ ) * tsA4x3_;
      tda4_   = ( tda4_ + tsa4_ ) * tsA4x3_;

      checkResult( tda4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "TDVecTSMatMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = ( tda4_ + tsa4_ ) * tsA4x3_;
      tsa4_   = ( tda4_ + tsa4_ ) * tsA4x3_;

      checkResult( tsa4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "TDVecTSMatMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ += tdc3_ * tsB3x3_;
      tdc3_   += tdc3_ * tsB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "TDVecTSMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ += ( tdc3_ + tsb3_ ) * tsB3x3_;
      tdc3_   += ( tdc3_ + tsb3_ ) * tsB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "TDVecTSMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ += ( tdc3_ + tsb3_ ) * tsB3x3_;
      tsb3_   += ( tdc3_ + tsb3_ ) * tsB3x3_;

      checkResult( tsb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "TDVecTSMatMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ -= tdc3_ * tsB3x3_;
      tdc3_   -= tdc3_ * tsB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "TDVecTSMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ -= ( tdc3_ + tsb3_ ) * tsB3x3_;
      tdc3_   -= ( tdc3_ + tsb3_ ) * tsB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "TDVecTSMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ -= ( tdc3_ + tsb3_ ) * tsB3x3_;
      tsb3_   -= ( tdc3_ + tsb3_ ) * tsB3x3_;

      checkResult( tsb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "TDVecTSMatMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ *= tdc3_ * tsB3x3_;
      tdc3_   *= tdc3_ * tsB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "TDVecTSMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ *= ( tdc3_ + tsb3_ ) * tsB3x3_;
      tdc3_   *= ( tdc3_ + tsb3_ ) * tsB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "TDVecTSMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ *= ( tdc3_ + tsb3_ ) * tsB3x3_;
      tsb3_   *= ( tdc3_ + tsb3_ ) * tsB3x3_;

      checkResult( tsb3_, result_ );
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
   // Initialization of the dense matrices
   //=====================================================================================

   // Initializing the first row-major dense matrix
   sA4x3_(0,0) = -1;
   sA4x3_(0,1) =  0;
   sA4x3_(0,2) = -2;
   sA4x3_(1,0) =  0;
   sA4x3_(1,1) =  2;
   sA4x3_(1,2) = -3;
   sA4x3_(2,0) =  0;
   sA4x3_(2,1) =  1;
   sA4x3_(2,2) =  2;
   sA4x3_(3,0) =  1;
   sA4x3_(3,1) =  0;
   sA4x3_(3,2) = -2;

   // Initializing the second row-major dense matrix
   sB3x3_(0,0) =  0;
   sB3x3_(0,0) = -1;
   sB3x3_(0,2) =  0;
   sB3x3_(1,1) =  1;
   sB3x3_(1,1) = -2;
   sB3x3_(1,1) =  2;
   sB3x3_(2,0) =  0;
   sB3x3_(2,0) =  0;
   sB3x3_(2,2) = -3;

   // Initializing the first column-major dense matrix
   tsA4x3_(0,0) = -1;
   tsA4x3_(0,1) =  0;
   tsA4x3_(0,2) = -2;
   tsA4x3_(1,0) =  0;
   tsA4x3_(1,1) =  2;
   tsA4x3_(1,2) = -3;
   tsA4x3_(2,0) =  0;
   tsA4x3_(2,1) =  1;
   tsA4x3_(2,2) =  2;
   tsA4x3_(3,0) =  1;
   tsA4x3_(3,1) =  0;
   tsA4x3_(3,2) = -2;

   // Initializing the second column-major dense matrix
   tsB3x3_(0,0) =  0;
   tsB3x3_(0,0) = -1;
   tsB3x3_(0,2) =  0;
   tsB3x3_(1,1) =  1;
   tsB3x3_(1,1) = -2;
   tsB3x3_(1,1) =  2;
   tsB3x3_(2,0) =  0;
   tsB3x3_(2,0) =  0;
   tsB3x3_(2,2) = -3;


   //=====================================================================================
   // Initialization of the dense vectors
   //=====================================================================================

   // Initializing the first dense row vector
   tda4_.resize( 4UL, false );
   tda4_[0] = -1;
   tda4_[1] =  0;
   tda4_[2] = -3;
   tda4_[3] =  2;

   // Initializing the second dense row vector
   tdb4_.resize( 4UL, false );
   tdb4_[0] =  0;
   tdb4_[1] =  1;
   tdb4_[2] =  2;
   tdb4_[3] = -1;

   // Initializing the third dense row vector
   tdc3_.resize( 3UL, false );
   tdc3_[0] = 1;
   tdc3_[1] = 2;
   tdc3_[2] = 3;

   // Initializing the fourth dense row vector
   tdd3_.resize( 3UL, false );
   tdd3_[0] = 0;
   tdd3_[1] = 2;
   tdd3_[2] = 1;


   //=====================================================================================
   // Initialization of the sparse vectors
   //=====================================================================================

   // Initializing the first sparse row vector
   tsa4_.resize( 4UL, false );
   tsa4_.reset();
   tsa4_[0] = -1;
   tsa4_[2] = -3;
   tsa4_[2] =  2;

   // Initializing the second sparse row vector
   tsb3_.resize( 3UL, false );
   tsb3_.reset();
   tsb3_[1] = 2;
   tsb3_[2] = 1;
}
//*************************************************************************************************

} // namespace tdvecsmatmult

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
      RUN_TDVECSMATMULT_ALIASING_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aliasing test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
