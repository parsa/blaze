//=================================================================================================
/*!
//  \file src/mathtest/rows/DenseGeneralTest1.cpp
//  \brief Source file for the Rows dense general test (part 1)
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <memory>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blazetest/mathtest/rows/DenseGeneralTest.h>


namespace blazetest {

namespace mathtest {

namespace rows {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Rows dense general test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseGeneralTest::DenseGeneralTest()
   : mat_ ( 5UL, 4UL )
   , tmat_( 5UL, 4UL )
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testSchurAssign();
   testMultAssign();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the Rows constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Rows specialization. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseGeneralTest::testConstructors()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Rows constructor";

      initialize();

      // Setup of empty row selection
      {
         std::vector<size_t> indices;
         RT rs = blaze::rows( mat_, indices.data(), 0UL );

         if( rs.rows() != 0UL || rs.columns() != mat_.columns() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty row selection failed\n"
                << " Details:\n"
                << "   Result:\n" << rs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of random in-bounds row selection
      for( size_t rep=0UL; rep<100UL; ++rep )
      {
         blaze::DynamicVector<size_t> indices( blaze::rand<size_t>( 1UL, 20UL ) );
         randomize( indices, 0UL, mat_.rows()-1UL );
         RT rs = blaze::rows( mat_, indices.data(), indices.size() );

         for( size_t i=0UL; i<rs.rows(); ++i ) {
            for( size_t j=0UL; j<rs.columns(); ++j ) {
               if( rs(i,j) != mat_(indices[i],j) ) {
                  std::ostringstream oss;
                  oss << " Test: " << test_ << "\n"
                      << " Error: Setup of row selection failed\n"
                      << " Details:\n"
                      << "   Indices:\n" << indices << "\n"
                      << "   Row selection:\n" << rs << "\n"
                      << "   Matrix:\n" << mat_ << "\n";
                  throw std::runtime_error( oss.str() );
               }
            }
         }
      }

      // Trying to setup an out-of-bounds row selection
      try {
         RT rs = blaze::rows( mat_, { 5 } );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds row selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Rows constructor";

      initialize();

      // Setup of empty row selection
      {
         std::vector<size_t> indices;
         ORT rs = blaze::rows( tmat_, indices.data(), 0UL );

         if( rs.rows() != 0UL || rs.columns() != tmat_.columns() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty row selection failed\n"
                << " Details:\n"
                << "   Result:\n" << rs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of random in-bounds row selection
      for( size_t rep=0UL; rep<100UL; ++rep )
      {
         blaze::DynamicVector<size_t> indices( blaze::rand<size_t>( 1UL, 20UL ) );
         randomize( indices, 0UL, tmat_.rows()-1UL );
         ORT rs = blaze::rows( tmat_, indices.data(), indices.size() );

         for( size_t i=0UL; i<rs.rows(); ++i ) {
            for( size_t j=0UL; j<rs.columns(); ++j ) {
               if( rs(i,j) != tmat_(indices[i],j) ) {
                  std::ostringstream oss;
                  oss << " Test: " << test_ << "\n"
                      << " Error: Setup of row selection failed\n"
                      << " Details:\n"
                      << "   Indices:\n" << indices << "\n"
                      << "   Row selection:\n" << rs << "\n"
                      << "   Matrix:\n" << tmat_ << "\n";
                  throw std::runtime_error( oss.str() );
               }
            }
         }
      }

      // Trying to setup an out-of-bounds row selection
      try {
         ORT rs = blaze::rows( tmat_, { 5 } );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds row selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Rows assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Rows specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseGeneralTest::testAssignment()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major Rows homogeneous assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );
      rs = 12;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 14UL );

      if( rs(0,0) != 12 || rs(0,1) != 12 || rs(0,2) != 12 || rs(0,3) != 12 ||
          rs(1,0) != 12 || rs(1,1) != 12 || rs(1,2) != 12 || rs(1,3) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 12 12 12 12 )\n( 12 12 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) != 12 || mat_(1,1) != 12 || mat_(1,2) != 12 || mat_(1,3) != 12 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 12 || mat_(3,1) != 12 || mat_(3,2) != 12 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 12 12 12 12 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 12 12 12 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major Rows list assignment (complete list)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );
      rs = { { 11, 0, 0, 12 }, { 0, 13, 14, 0 } };

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Rows list assignment (incomplete list)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );
      rs = { { 11, 0, 0, 12 }, { 0, 13, 14 } };

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major Rows copy assignment (no aliasing)";

      initialize();

      MT mat{ {  0,  0,  0,  0 },
              { 11,  0, 12,  0 },
              {  0,  0,  0,  0 },
              { 13, 14, 15, 16 },
              {  0,  0,  0,  0 } };

      RT rs = blaze::rows( mat, { 3UL, 1UL } );
      rs = blaze::rows( mat_, { 3UL, 1UL } );

      checkRows    ( rs , 2UL );
      checkColumns ( rs , 4UL );
      checkNonZeros( rs , 4UL );
      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( rs(0,0) != 0 || rs(0,1) != 4 || rs(0,2) != 5 || rs(0,3) != -6 ||
          rs(1,0) != 0 || rs(1,1) != 1 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(  0  4  5 -6 )\n(  0  1  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) != 1 || mat(1,2) != 0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 4 || mat(3,2) != 5 || mat(3,3) != -6 ||
          mat(4,0) != 0 || mat(4,1) != 0 || mat(4,2) != 0 || mat(4,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  0  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Rows copy assignment (aliasing)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 4UL } );
      rs = blaze::rows( mat_, { 2UL, 3UL } );

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 5UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 8UL );

      if( rs(0,0) != -2 || rs(0,1) != 0 || rs(0,2) != -3 || rs(0,3) !=  0 ||
          rs(1,0) !=  0 || rs(1,1) != 4 || rs(1,2) !=  5 || rs(1,3) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -2  0 -3  0 )\n(  0  4  5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != -2 || mat_(3,1) !=  0 || mat_(3,2) != -3 || mat_(3,3) !=  0 ||
          mat_(4,0) !=  0 || mat_(4,1) !=  4 || mat_(4,2) !=  5 || mat_(4,3) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { 11,  0,  0, 12 },
                                                      {  0, 13, 14,  0 } };

      rs = mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs = mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs = mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { 11,  0,  0, 12 },
                                                         {  0, 13, 14,  0 } };

      rs = mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs = mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs = mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ { 11,  0,  0, 12 },
                                                       {  0, 13, 14,  0 } };

      rs = mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ { 11,  0,  0, 12 },
                                                          {  0, 13, 14,  0 } };

      rs = mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  4UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 13 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major Rows homogeneous assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );
      rs = 12;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 14UL );

      if( rs(0,0) != 12 || rs(0,1) != 12 || rs(0,2) != 12 || rs(0,3) != 12 ||
          rs(1,0) != 12 || rs(1,1) != 12 || rs(1,2) != 12 || rs(1,3) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 12 12 12 12 )\n( 12 12 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 12 || tmat_(1,1) != 12 || tmat_(1,2) != 12 || tmat_(1,3) != 12 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 12 || tmat_(3,1) != 12 || tmat_(3,2) != 12 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 12 12 12 12 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 12 12 12 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major Rows list assignment (complete list)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );
      rs = { { 11, 0, 0, 12 }, { 0, 13, 14, 0 } };

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Rows list assignment (incomplete list)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );
      rs = { { 11, 0, 0, 12 }, { 0, 13, 14 } };

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major Rows copy assignment (no aliasing)";

      initialize();

      OMT mat{ {  0,  0,  0,  0 },
               { 11,  0, 12,  0 },
               {  0,  0,  0,  0 },
               { 13, 14, 15, 16 },
               {  0,  0,  0,  0 } };

      ORT rs = blaze::rows( mat, { 3UL, 1UL } );
      rs = blaze::rows( tmat_, { 3UL, 1UL } );

      checkRows    ( rs , 2UL );
      checkColumns ( rs , 4UL );
      checkNonZeros( rs , 4UL );
      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( rs(0,0) != 0 || rs(0,1) != 4 || rs(0,2) != 5 || rs(0,3) != -6 ||
          rs(1,0) != 0 || rs(1,1) != 1 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(  0  4  5 -6 )\n(  0  1  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) != 1 || mat(1,2) != 0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 4 || mat(3,2) != 5 || mat(3,3) != -6 ||
          mat(4,0) != 0 || mat(4,1) != 0 || mat(4,2) != 0 || mat(4,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  0  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Rows copy assignment (aliasing)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 4UL } );
      rs = blaze::rows( tmat_, { 2UL, 3UL } );

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 5UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 8UL );

      if( rs(0,0) != -2 || rs(0,1) != 0 || rs(0,2) != -3 || rs(0,3) !=  0 ||
          rs(1,0) !=  0 || rs(1,1) != 4 || rs(1,2) !=  5 || rs(1,3) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -2  0 -3  0 )\n(  0  4  5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -2 || tmat_(3,1) !=  0 || tmat_(3,2) != -3 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  4 || tmat_(4,2) !=  5 || tmat_(4,3) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0  4  5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix assignment (mixed type)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { 11,  0,  0, 12 },
                                                      {  0, 13, 14,  0 } };

      rs = mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix assignment (aligned/padded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs = mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs = mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (mixed type)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { 11,  0,  0, 12 },
                                                         {  0, 13, 14,  0 } };

      rs = mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (aligned/padded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs = mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs = mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ { 11,  0,  0, 12 },
                                                       {  0, 13, 14,  0 } };

      rs = mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ { 11,  0,  0, 12 },
                                                          {  0, 13, 14,  0 } };

      rs = mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  4UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( rs(0,0) != 11 || rs(0,1) !=  0 || rs(0,2) !=  0 || rs(0,3) != 12 ||
          rs(1,0) !=  0 || rs(1,1) != 13 || rs(1,2) != 14 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  0  0 12 )\n(  0 13 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 13 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 13 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  0  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Rows addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Rows specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseGeneralTest::testAddAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Rows addition assignment
   //=====================================================================================

   {
      test_ = "Row-major Rows addition assignment (no aliasing)";

      initialize();

      MT mat{ {  0,  0,  0,  0 },
              { 11,  0, 12,  0 },
              {  0,  0,  0,  0 },
              { 13, 14, 15, 16 },
              {  0,  0,  0,  0 } };

      RT rs = blaze::rows( mat, { 3UL, 1UL } );
      rs += blaze::rows( mat_, { 3UL, 1UL } );

      checkRows    ( rs , 2UL );
      checkColumns ( rs , 4UL );
      checkNonZeros( rs , 7UL );
      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 7UL );

      if( rs(0,0) != 13 || rs(0,1) != 18 || rs(0,2) != 20 || rs(0,3) != 10 ||
          rs(1,0) != 11 || rs(1,1) !=  1 || rs(1,2) != 12 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 13 18 20 10 )\n( 11  1 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=  0 || mat(0,1) !=  0 || mat(0,2) !=  0 || mat(0,3) !=  0 ||
          mat(1,0) != 11 || mat(1,1) !=  1 || mat(1,2) != 12 || mat(1,3) !=  0 ||
          mat(2,0) !=  0 || mat(2,1) !=  0 || mat(2,2) !=  0 || mat(2,3) !=  0 ||
          mat(3,0) != 13 || mat(3,1) != 18 || mat(3,2) != 20 || mat(3,3) != 10 ||
          mat(4,0) !=  0 || mat(4,1) !=  0 || mat(4,2) !=  0 || mat(4,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 11  1 12  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "( 13 18 20 10 )\n"
                                     "(  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Rows addition assignment (aliasing)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 4UL } );
      rs += blaze::rows( mat_, { 2UL, 3UL } );

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 8UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 11UL );

      if( rs(0,0) != -2 || rs(0,1) !=  4 || rs(0,2) !=  2 || rs(0,3) != -6 ||
          rs(1,0) !=  7 || rs(1,1) != -4 || rs(1,2) != 14 || rs(1,3) !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -2  4  2 -6 )\n(  7 -4 14  4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != -2 || mat_(3,1) !=  4 || mat_(3,2) !=  2 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -4 || mat_(4,2) != 14 || mat_(4,3) !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( -2  4  2 -6 )\n"
                                     "(  7 -4 14  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix addition assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { 11,  0,  0, 12 },
                                                      {  0, 13, 14,  0 } };

      rs += mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 14 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) !=  6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix addition assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs += mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 14 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) !=  6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs += mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 14 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) !=  6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { 11,  0,  0, 12 },
                                                         {  0, 13, 14,  0 } };

      rs += mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 14 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) !=  6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs += mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 14 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) !=  6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs += mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 14 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) !=  6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix addition assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ { 11,  0,  0, 12 },
                                                       {  0, 13, 14,  0 } };

      rs += mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 14 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) !=  6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix addition assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ { 11,  0,  0, 12 },
                                                          {  0, 13, 14,  0 } };

      rs += mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 14 || mat_(1,2) != 14 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) != 11 || mat_(3,1) !=  4 || mat_(3,2) !=  5 || mat_(3,3) !=  6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Rows addition assignment
   //=====================================================================================

   {
      test_ = "Column-major Rows addition assignment (no aliasing)";

      initialize();

      OMT mat{ {  0,  0,  0,  0 },
               { 11,  0, 12,  0 },
               {  0,  0,  0,  0 },
               { 13, 14, 15, 16 },
               {  0,  0,  0,  0 } };

      ORT rs = blaze::rows( mat, { 3UL, 1UL } );
      rs += blaze::rows( tmat_, { 3UL, 1UL } );

      checkRows    ( rs , 2UL );
      checkColumns ( rs , 4UL );
      checkNonZeros( rs , 7UL );
      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 7UL );

      if( rs(0,0) != 13 || rs(0,1) != 18 || rs(0,2) != 20 || rs(0,3) != 10 ||
          rs(1,0) != 11 || rs(1,1) !=  1 || rs(1,2) != 12 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 13 18 20 10 )\n( 11  1 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=  0 || mat(0,1) !=  0 || mat(0,2) !=  0 || mat(0,3) !=  0 ||
          mat(1,0) != 11 || mat(1,1) !=  1 || mat(1,2) != 12 || mat(1,3) !=  0 ||
          mat(2,0) !=  0 || mat(2,1) !=  0 || mat(2,2) !=  0 || mat(2,3) !=  0 ||
          mat(3,0) != 13 || mat(3,1) != 18 || mat(3,2) != 20 || mat(3,3) != 10 ||
          mat(4,0) !=  0 || mat(4,1) !=  0 || mat(4,2) !=  0 || mat(4,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 11  1 12  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "( 13 18 20 10 )\n"
                                     "(  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Rows addition assignment (aliasing)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 4UL } );
      rs += blaze::rows( tmat_, { 2UL, 3UL } );

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 8UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 11UL );

      if( rs(0,0) != -2 || rs(0,1) !=  4 || rs(0,2) !=  2 || rs(0,3) != -6 ||
          rs(1,0) !=  7 || rs(1,1) != -4 || rs(1,2) != 14 || rs(1,3) !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -2  4  2 -6 )\n(  7 -4 14  4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -2 || tmat_(3,1) !=  4 || tmat_(3,2) !=  2 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -4 || tmat_(4,2) != 14 || tmat_(4,3) !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  1  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( -2  4  2 -6 )\n"
                                     "(  7 -4 14  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix addition assignment (mixed type)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { 11,  0,  0, 12 },
                                                      {  0, 13, 14,  0 } };

      rs += mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  6UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 14 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix addition assignment (aligned/padded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs += mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  6UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 14 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs += mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  6UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 14 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (mixed type)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { 11,  0,  0, 12 },
                                                         {  0, 13, 14,  0 } };

      rs += mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  6UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 14 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (aligned/padded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs += mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  6UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 14 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs += mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  6UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 14 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix addition assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ { 11,  0,  0, 12 },
                                                       {  0, 13, 14,  0 } };

      rs += mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  6UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 14 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix addition assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ { 11,  0,  0, 12 },
                                                          {  0, 13, 14,  0 } };

      rs += mat;

      checkRows    ( rs   ,  2UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  6UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) != 11 || rs(0,1) !=  4 || rs(0,2) !=  5 || rs(0,3) != 6 ||
          rs(1,0) !=  0 || rs(1,1) != 14 || rs(1,2) != 14 || rs(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 11  4 17 -6 )\n(  0 14 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 14 || tmat_(1,2) != 14 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 11 || tmat_(3,1) !=  4 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0 14 14  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "( 11  4  5  6 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Rows subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Rows
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseGeneralTest::testSubAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Rows subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major Rows subtraction assignment (no aliasing)";

      initialize();

      MT mat{ {  0,  0,  0,  0 },
              { 11,  0, 12,  0 },
              {  0,  0,  0,  0 },
              { 13, 14, 15, 16 },
              {  0,  0,  0,  0 } };

      RT rs = blaze::rows( mat, { 3UL, 1UL } );
      rs -= blaze::rows( mat_, { 3UL, 1UL } );

      checkRows    ( rs , 2UL );
      checkColumns ( rs , 4UL );
      checkNonZeros( rs , 7UL );
      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 7UL );

      if( rs(0,0) != 13 || rs(0,1) != 10 || rs(0,2) != 10 || rs(0,3) != 22 ||
          rs(1,0) != 11 || rs(1,1) != -1 || rs(1,2) != 12 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 13 10 10 22 )\n( 11 -1 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=  0 || mat(0,1) !=  0 || mat(0,2) !=  0 || mat(0,3) !=  0 ||
          mat(1,0) != 11 || mat(1,1) != -1 || mat(1,2) != 12 || mat(1,3) !=  0 ||
          mat(2,0) !=  0 || mat(2,1) !=  0 || mat(2,2) !=  0 || mat(2,3) !=  0 ||
          mat(3,0) != 13 || mat(3,1) != 10 || mat(3,2) != 10 || mat(3,3) != 22 ||
          mat(4,0) !=  0 || mat(4,1) !=  0 || mat(4,2) !=  0 || mat(4,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 11 -1 12  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "( 13 10 10 22 )\n"
                                     "(  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Rows subtraction assignment (aliasing)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 4UL } );
      rs -= blaze::rows( mat_, { 2UL, 3UL } );

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 8UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 11UL );

      if( rs(0,0) != 2 || rs(0,1) !=   4 || rs(0,2) != 8 || rs(0,3) != -6 ||
          rs(1,0) != 7 || rs(1,1) != -12 || rs(1,2) != 4 || rs(1,3) != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 2   4  8 -6 )\n( 7 -12  4 16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=   0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=   1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=   0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  2 || mat_(3,1) !=   4 || mat_(3,2) !=  8 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -12 || mat_(4,2) !=  4 || mat_(4,3) != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0   0  0  0 )\n"
                                     "(  0   1  0  0 )\n"
                                     "( -2   0 -3  0 )\n"
                                     "(  2   4  8 -6 )\n"
                                     "(  7 -12  4 16 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { 11,  0,  0, 12 },
                                                      {  0, 13, 14,  0 } };

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { 11,  0,  0, 12 },
                                                         {  0, 13, 14,  0 } };

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix subtraction assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ { 11,  0,  0, 12 },
                                                       {  0, 13, 14,  0 } };

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix subtraction assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ { 11,  0,  0, 12 },
                                                          {  0, 13, 14,  0 } };

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Rows subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major Rows subtraction assignment (no aliasing)";

      initialize();

      MT mat{ {  0,  0,  0,  0 },
              { 11,  0, 12,  0 },
              {  0,  0,  0,  0 },
              { 13, 14, 15, 16 },
              {  0,  0,  0,  0 } };

      RT rs = blaze::rows( mat, { 3UL, 1UL } );
      rs -= blaze::rows( mat_, { 3UL, 1UL } );

      checkRows    ( rs , 2UL );
      checkColumns ( rs , 4UL );
      checkNonZeros( rs , 7UL );
      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 7UL );

      if( rs(0,0) != 13 || rs(0,1) != 10 || rs(0,2) != 10 || rs(0,3) != 22 ||
          rs(1,0) != 11 || rs(1,1) != -1 || rs(1,2) != 12 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 13 10 10 22 )\n( 11 -1 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=  0 || mat(0,1) !=  0 || mat(0,2) !=  0 || mat(0,3) !=  0 ||
          mat(1,0) != 11 || mat(1,1) != -1 || mat(1,2) != 12 || mat(1,3) !=  0 ||
          mat(2,0) !=  0 || mat(2,1) !=  0 || mat(2,2) !=  0 || mat(2,3) !=  0 ||
          mat(3,0) != 13 || mat(3,1) != 10 || mat(3,2) != 10 || mat(3,3) != 22 ||
          mat(4,0) !=  0 || mat(4,1) !=  0 || mat(4,2) !=  0 || mat(4,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "( 11 -1 12  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "( 13 10 10 22 )\n"
                                     "(  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Rows subtraction assignment (aliasing)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 4UL } );
      rs -= blaze::rows( mat_, { 2UL, 3UL } );

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 8UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 11UL );

      if( rs(0,0) != 2 || rs(0,1) !=   4 || rs(0,2) != 8 || rs(0,3) != -6 ||
          rs(1,0) != 7 || rs(1,1) != -12 || rs(1,2) != 4 || rs(1,3) != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 2   4  8 -6 )\n( 7 -12  4 16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=   0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=   1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=   0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  2 || mat_(3,1) !=   4 || mat_(3,2) !=  8 || mat_(3,3) != -6 ||
          mat_(4,0) !=  7 || mat_(4,1) != -12 || mat_(4,2) !=  4 || mat_(4,3) != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0   0  0  0 )\n"
                                     "(  0   1  0  0 )\n"
                                     "( -2   0 -3  0 )\n"
                                     "(  2   4  8 -6 )\n"
                                     "(  7 -12  4 16 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { 11,  0,  0, 12 },
                                                      {  0, 13, 14,  0 } };

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { 11,  0,  0, 12 },
                                                         {  0, 13, 14,  0 } };

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,0) = 11;
      mat(0,3) = 12;
      mat(1,1) = 13;
      mat(1,2) = 14;

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix subtraction assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ { 11,  0,  0, 12 },
                                                       {  0, 13, 14,  0 } };

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix subtraction assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ { 11,  0,  0, 12 },
                                                          {  0, 13, 14,  0 } };

      rs -= mat;

      checkRows    ( rs  ,  2UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  6UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) != -11 || rs(0,1) !=   4 || rs(0,2) !=   5 || rs(0,3) != -18 ||
          rs(1,0) !=   0 || rs(1,1) != -12 || rs(1,2) != -14 || rs(1,3) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( -11   4   5 -18 )\n(   0 -12 -14   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=   0 || mat_(1,1) != -12 || mat_(1,2) != -14 || mat_(1,3) !=   0 ||
          mat_(2,0) !=  -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) != -11 || mat_(3,1) !=   4 || mat_(3,2) !=   5 || mat_(3,3) != -18 ||
          mat_(4,0) !=   7 || mat_(4,1) !=  -8 || mat_(4,2) !=   9 || mat_(4,3) !=  10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(   0 -12 -14   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "( -11   4   5 -18 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Rows Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the Rows
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseGeneralTest::testSchurAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Rows Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major Rows Schur product assignment (no aliasing)";

      initialize();

      MT mat{ { 0, 0, 0, 0 },
              { 1, 2, 3, 0 },
              { 0, 0, 0, 0 },
              { 4, 3, 2, 1 },
              { 0, 0, 0, 0 } };

      RT rs = blaze::rows( mat, { 3UL, 1UL } );
      rs %= blaze::rows( mat_, { 3UL, 1UL } );

      checkRows    ( rs , 2UL );
      checkColumns ( rs , 4UL );
      checkNonZeros( rs , 4UL );
      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( rs(0,0) != 0 || rs(0,1) != 12 || rs(0,2) != 10 || rs(0,3) != -6 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) !=  0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 12 10 -6 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  0 || mat(0,2) !=  0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) !=  2 || mat(1,2) !=  0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) !=  0 || mat(2,2) !=  0 || mat(2,3) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 12 || mat(3,2) != 10 || mat(3,3) != -6 ||
          mat(4,0) != 0 || mat(4,1) !=  0 || mat(4,2) !=  0 || mat(4,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "(  0 12 10 -6 )\n"
                                     "(  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Rows Schur product assignment (aliasing)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 4UL } );
      rs %= blaze::rows( mat_, { 2UL, 3UL } );

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 4UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 7UL );

      if( rs(0,0) != 0 || rs(0,1) !=   0 || rs(0,2) != -15 || rs(0,3) !=   0 ||
          rs(1,0) != 0 || rs(1,1) != -32 || rs(1,2) !=  45 || rs(1,3) != -60 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0   0 -15   0 )\n( 0 -32  45 -60 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=   0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=   1 || mat_(1,2) !=   0 || mat_(1,3) !=   0 ||
          mat_(2,0) != -2 || mat_(2,1) !=   0 || mat_(2,2) !=  -3 || mat_(2,3) !=   0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=   0 || mat_(3,2) != -15 || mat_(3,3) !=   0 ||
          mat_(4,0) !=  0 || mat_(4,1) != -32 || mat_(4,2) !=  45 || mat_(4,3) != -60 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0   0   0   0 )\n"
                                     "(  0   1   0   0 )\n"
                                     "( -2   0  -3   0 )\n"
                                     "(  0   0 -15   0 )\n"
                                     "(  0 -32  45 -60 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix Schur product assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { 0, -1, 0, -2 },
                                                      { 0,  2, 1,  0 } };

      rs %= mat;

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != -4 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix Schur product assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,1) = -1;
      mat(0,3) = -2;
      mat(1,1) =  2;
      mat(1,2) =  1;

      rs %= mat;

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != -4 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix Schur product assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,1) = -1;
      mat(0,3) = -2;
      mat(1,1) =  2;
      mat(1,2) =  1;

      rs %= mat;

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != -4 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix Schur product assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { 0, -1, 0, -2 },
                                                         { 0,  2, 1,  0 } };

      rs %= mat;

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != -4 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix Schur product assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,1) = -1;
      mat(0,3) = -2;
      mat(1,1) =  2;
      mat(1,2) =  1;

      rs %= mat;

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != -4 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix Schur product assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,1) = -1;
      mat(0,3) = -2;
      mat(1,1) =  2;
      mat(1,2) =  1;

      rs %= mat;

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != -4 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix Schur product assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ { 0, -1, 0, -2 },
                                                       { 0,  2, 1,  0 } };

      rs %= mat;

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != -4 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix Schur product assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ { 0, -1, 0, -2 },
                                                          { 0,  2, 1,  0 } };

      rs %= mat;

      checkRows    ( rs  , 2UL );
      checkColumns ( rs  , 4UL );
      checkNonZeros( rs  , 3UL );
      checkRows    ( mat_, 5UL );
      checkColumns ( mat_, 4UL );
      checkNonZeros( mat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=  0 || mat_(0,1) !=  0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) !=  2 || mat_(1,2) !=  0 || mat_(1,3) !=  0 ||
          mat_(2,0) != -2 || mat_(2,1) !=  0 || mat_(2,2) != -3 || mat_(2,3) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != -4 || mat_(3,2) !=  0 || mat_(3,3) != 12 ||
          mat_(4,0) !=  7 || mat_(4,1) != -8 || mat_(4,2) !=  9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Rows Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major Rows Schur product assignment (no aliasing)";

      initialize();

      OMT mat{ { 0, 0, 0, 0 },
               { 1, 2, 3, 0 },
               { 0, 0, 0, 0 },
               { 4, 3, 2, 1 },
               { 0, 0, 0, 0 } };

      ORT rs = blaze::rows( mat, { 3UL, 1UL } );
      rs %= blaze::rows( tmat_, { 3UL, 1UL } );

      checkRows    ( rs , 2UL );
      checkColumns ( rs , 4UL );
      checkNonZeros( rs , 4UL );
      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 4UL );

      if( rs(0,0) != 0 || rs(0,1) != 12 || rs(0,2) != 10 || rs(0,3) != -6 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) !=  0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 12 10 -6 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  0 || mat(0,2) !=  0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) !=  2 || mat(1,2) !=  0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) !=  0 || mat(2,2) !=  0 || mat(2,3) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 12 || mat(3,2) != 10 || mat(3,3) != -6 ||
          mat(4,0) != 0 || mat(4,1) !=  0 || mat(4,2) !=  0 || mat(4,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "(  0 12 10 -6 )\n"
                                     "(  0  0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Rows Schur product assignment (aliasing)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 4UL } );
      rs %= blaze::rows( tmat_, { 2UL, 3UL } );

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 4UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 7UL );

      if( rs(0,0) != 0 || rs(0,1) !=   0 || rs(0,2) != -15 || rs(0,3) !=   0 ||
          rs(1,0) != 0 || rs(1,1) != -32 || rs(1,2) !=  45 || rs(1,3) != -60 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0   0 -15   0 )\n( 0 -32  45 -60 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=   0 || tmat_(0,2) !=   0 || tmat_(0,3) !=   0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=   1 || tmat_(1,2) !=   0 || tmat_(1,3) !=   0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=   0 || tmat_(2,2) !=  -3 || tmat_(2,3) !=   0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) !=   0 || tmat_(3,2) != -15 || tmat_(3,3) !=   0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) != -32 || tmat_(4,2) !=  45 || tmat_(4,3) != -60 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0   0   0   0 )\n"
                                     "(  0   1   0   0 )\n"
                                     "( -2   0  -3   0 )\n"
                                     "(  0   0 -15   0 )\n"
                                     "(  0 -32  45 -60 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix Schur product assignment (mixed type)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ { 0, -1, 0, -2 },
                                                      { 0,  2, 1,  0 } };

      rs %= mat;

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  2 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -4 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix Schur product assignment (aligned/padded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,1) = -1;
      mat(0,3) = -2;
      mat(1,1) =  2;
      mat(1,2) =  1;

      rs %= mat;

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  2 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -4 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix Schur product assignment (unaligned/unpadded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,1) = -1;
      mat(0,3) = -2;
      mat(1,1) =  2;
      mat(1,2) =  1;

      rs %= mat;

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  2 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -4 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix Schur product assignment (mixed type)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ { 0, -1, 0, -2 },
                                                         { 0,  2, 1,  0 } };

      rs %= mat;

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  2 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -4 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix Schur product assignment (aligned/padded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 2UL, 4UL, 16UL );
      mat = 0;
      mat(0,1) = -1;
      mat(0,3) = -2;
      mat(1,1) =  2;
      mat(1,2) =  1;

      rs %= mat;

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  2 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -4 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix Schur product assignment (unaligned/unpadded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[9UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 2UL, 4UL );
      mat = 0;
      mat(0,1) = -1;
      mat(0,3) = -2;
      mat(1,1) =  2;
      mat(1,2) =  1;

      rs %= mat;

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  2 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -4 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix Schur product assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ { 0, -1, 0, -2 },
                                                       { 0,  2, 1,  0 } };

      rs %= mat;

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  2 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -4 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix Schur product assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ { 0, -1, 0, -2 },
                                                          { 0,  2, 1,  0 } };

      rs %= mat;

      checkRows    ( rs   , 2UL );
      checkColumns ( rs   , 4UL );
      checkNonZeros( rs   , 3UL );
      checkRows    ( tmat_, 5UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( rs(0,0) != 0 || rs(0,1) != -4 || rs(0,2) != 0 || rs(0,3) != 12 ||
          rs(1,0) != 0 || rs(1,1) !=  2 || rs(1,2) != 0 || rs(1,3) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n( 0 -4  0 12 )\n( 0  2  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=  0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  2 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != -2 || tmat_(2,1) !=  0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -4 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 ||
          tmat_(4,0) !=  7 || tmat_(4,1) != -8 || tmat_(4,2) !=  9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(  0  0  0  0 )\n"
                                     "(  0  2  0  0 )\n"
                                     "( -2  0 -3  0 )\n"
                                     "(  0 -4  0 12 )\n"
                                     "(  7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Rows multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Rows
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseGeneralTest::testMultAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Rows multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major Rows multiplication assignment (no aliasing)";

      initialize();

      MT mat{ {  0,  0,  0,  0 },
              {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              {  0,  4,  5, -6 },
              {  7, -8,  9, 10 } };

      RT rs = blaze::rows( mat, { 2UL, 0UL, 3UL, 1UL } );
      rs *= blaze::rows( mat_, { 1UL, 2UL, 2UL, 1UL } );

      checkRows    ( rs ,  4UL );
      checkColumns ( rs ,  4UL );
      checkNonZeros( rs ,  8UL );
      checkRows    ( mat,  5UL );
      checkColumns ( mat,  4UL );
      checkNonZeros( mat, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=   0 || mat(0,1) !=  0 || mat(0,2) !=   0 || mat(0,3) !=  0 ||
          mat(1,0) !=  -2 || mat(1,1) !=  0 || mat(1,2) !=  -3 || mat(1,3) !=  0 ||
          mat(2,0) !=   6 || mat(2,1) != -2 || mat(2,2) !=   9 || mat(2,3) !=  0 ||
          mat(3,0) != -18 || mat(3,1) != -6 || mat(3,2) != -27 || mat(3,3) !=  0 ||
          mat(4,0) !=   7 || mat(4,1) != -8 || mat(4,2) !=   9 || mat(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Rows multiplication assignment (aliasing)";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );
      rs *= blaze::rows( mat_, { 1UL, 2UL, 2UL, 1UL } );

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {  0,  1,  0,  0 },
                                                      { -2,  0, -3,  0 },
                                                      { -2,  0, -3,  0 },
                                                      {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 4UL, 4UL, 16UL );
      mat = { {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              { -2,  0, -3,  0 },
              {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[17UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 4UL, 4UL );
      mat = { {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              { -2,  0, -3,  0 },
              {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (mixed type)";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {  0,  1,  0,  0 },
                                                         { -2,  0, -3,  0 },
                                                         { -2,  0, -3,  0 },
                                                         {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 4UL, 4UL, 16UL );
      mat = { {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              { -2,  0, -3,  0 },
              {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[17UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 4UL, 4UL );
      mat = { {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              { -2,  0, -3,  0 },
              {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix multiplication assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ {  0,  1,  0,  0 },
                                                       { -2,  0, -3,  0 },
                                                       { -2,  0, -3,  0 },
                                                       {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix multiplication assignment";

      initialize();

      RT rs = blaze::rows( mat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ {  0,  1,  0,  0 },
                                                          { -2,  0, -3,  0 },
                                                          { -2,  0, -3,  0 },
                                                          {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs  ,  4UL );
      checkColumns ( rs  ,  4UL );
      checkNonZeros( rs  ,  8UL );
      checkRows    ( mat_,  5UL );
      checkColumns ( mat_,  4UL );
      checkNonZeros( mat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   0 || mat_(0,1) !=  0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 ||
          mat_(1,0) !=  -2 || mat_(1,1) !=  0 || mat_(1,2) !=  -3 || mat_(1,3) !=  0 ||
          mat_(2,0) !=   6 || mat_(2,1) != -2 || mat_(2,2) !=   9 || mat_(2,3) !=  0 ||
          mat_(3,0) != -18 || mat_(3,1) != -6 || mat_(3,2) != -27 || mat_(3,3) !=  0 ||
          mat_(4,0) !=   7 || mat_(4,1) != -8 || mat_(4,2) !=   9 || mat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Rows multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major Rows multiplication assignment (no aliasing)";

      initialize();

      OMT mat{ {  0,  0,  0,  0 },
               {  0,  1,  0,  0 },
               { -2,  0, -3,  0 },
               {  0,  4,  5, -6 },
               {  7, -8,  9, 10 } };

      ORT rs = blaze::rows( mat, { 2UL, 0UL, 3UL, 1UL } );
      rs *= blaze::rows( tmat_, { 1UL, 2UL, 2UL, 1UL } );

      checkRows    ( rs ,  4UL );
      checkColumns ( rs ,  4UL );
      checkNonZeros( rs ,  8UL );
      checkRows    ( mat,  5UL );
      checkColumns ( mat,  4UL );
      checkNonZeros( mat, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=   0 || mat(0,1) !=  0 || mat(0,2) !=   0 || mat(0,3) !=  0 ||
          mat(1,0) !=  -2 || mat(1,1) !=  0 || mat(1,2) !=  -3 || mat(1,3) !=  0 ||
          mat(2,0) !=   6 || mat(2,1) != -2 || mat(2,2) !=   9 || mat(2,3) !=  0 ||
          mat(3,0) != -18 || mat(3,1) != -6 || mat(3,2) != -27 || mat(3,3) !=  0 ||
          mat(4,0) !=   7 || mat(4,1) != -8 || mat(4,2) !=   9 || mat(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Rows multiplication assignment (aliasing)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );
      rs *= blaze::rows( tmat_, { 1UL, 2UL, 2UL, 1UL } );

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (mixed type)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::DynamicMatrix<short,rowMajor> mat{ {  0,  1,  0,  0 },
                                                      { -2,  0, -3,  0 },
                                                      { -2,  0, -3,  0 },
                                                      {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 4UL, 4UL, 16UL );
      mat = { {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              { -2,  0, -3,  0 },
              {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[17UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 4UL, 4UL );
      mat = { {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              { -2,  0, -3,  0 },
              {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (mixed type)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::DynamicMatrix<short,columnMajor> mat{ {  0,  1,  0,  0 },
                                                         { -2,  0, -3,  0 },
                                                         { -2,  0, -3,  0 },
                                                         {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat( memory.get(), 4UL, 4UL, 16UL );
      mat = { {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              { -2,  0, -3,  0 },
              {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[17UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 4UL, 4UL );
      mat = { {  0,  1,  0,  0 },
              { -2,  0, -3,  0 },
              { -2,  0, -3,  0 },
              {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix multiplication assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::CompressedMatrix<int,rowMajor> mat{ {  0,  1,  0,  0 },
                                                       { -2,  0, -3,  0 },
                                                       { -2,  0, -3,  0 },
                                                       {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix multiplication assignment";

      initialize();

      ORT rs = blaze::rows( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::CompressedMatrix<int,columnMajor> mat{ {  0,  1,  0,  0 },
                                                          { -2,  0, -3,  0 },
                                                          { -2,  0, -3,  0 },
                                                          {  0,  1,  0,  0 } };

      rs *= mat;

      checkRows    ( rs   ,  4UL );
      checkColumns ( rs   ,  4UL );
      checkNonZeros( rs   ,  8UL );
      checkRows    ( tmat_,  5UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( rs(0,0) !=   6 || rs(0,1) != -2 || rs(0,2) !=   9 || rs(0,3) != 0 ||
          rs(1,0) !=   0 || rs(1,1) !=  0 || rs(1,2) !=   0 || rs(1,3) != 0 ||
          rs(2,0) != -18 || rs(2,1) != -6 || rs(2,2) != -27 || rs(2,3) != 0 ||
          rs(3,0) !=  -2 || rs(3,1) !=  0 || rs(3,2) !=  -3 || rs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << rs << "\n"
             << "   Expected result:\n(   6  -2   9  0 )\n"
                                     "(   0   0  -3  0 )\n"
                                     "( -18  -6 -27  0 )\n"
                                     "(  -2   0  -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  -2 || tmat_(1,1) !=  0 || tmat_(1,2) !=  -3 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=   6 || tmat_(2,1) != -2 || tmat_(2,2) !=   9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != -18 || tmat_(3,1) != -6 || tmat_(3,2) != -27 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=   7 || tmat_(4,1) != -8 || tmat_(4,2) !=   9 || tmat_(4,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   0   0   0   0 )\n"
                                     "(  -2   0  -3   0 )\n"
                                     "(   6  -2   9   0 )\n"
                                     "( -18  -6 -27   0 )\n"
                                     "(   7  -8   9  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member matrices to specific predetermined values.
*/
void DenseGeneralTest::initialize()
{
   // Initializing the row-major dynamic matrix
   mat_.reset();
   mat_(1,1) =  1;
   mat_(2,0) = -2;
   mat_(2,2) = -3;
   mat_(3,1) =  4;
   mat_(3,2) =  5;
   mat_(3,3) = -6;
   mat_(4,0) =  7;
   mat_(4,1) = -8;
   mat_(4,2) =  9;
   mat_(4,3) = 10;

   // Initializing the column-major dynamic matrix
   tmat_.reset();
   tmat_(1,1) =  1;
   tmat_(2,0) = -2;
   tmat_(2,2) = -3;
   tmat_(3,1) =  4;
   tmat_(3,2) =  5;
   tmat_(3,3) = -6;
   tmat_(4,0) =  7;
   tmat_(4,1) = -8;
   tmat_(4,2) =  9;
   tmat_(4,3) = 10;
}
//*************************************************************************************************

} // namespace rows

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
   std::cout << "   Running Rows dense general test (part 1)..." << std::endl;

   try
   {
      RUN_ROWS_DENSEGENERAL_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Rows dense general test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
