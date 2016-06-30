//=================================================================================================
/*!
//  \file src/mathtest/symmetricmatrix/DenseNonNumericTest.cpp
//  \brief Source file for the SymmetricMatrix dense non-numeric test
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
//  TO, PROCUREMENT OF SUSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
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
#include <blaze/math/DenseColumn.h>
#include <blaze/math/DenseRow.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/util/policies/ArrayDelete.h>
#include <blazetest/mathtest/symmetricmatrix/DenseNonNumericTest.h>


namespace blazetest {

namespace mathtest {

namespace symmetricmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SymmetricMatrix dense non-numeric test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseNonNumericTest::DenseNonNumericTest()
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
   testResize();
   testExtend();
   testReserve();
   testTranspose();
   testCTranspose();
   testSwap();
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
void DenseNonNumericTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Row-major SymmetricMatrix default constructor (StaticMatrix)";

      const blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> > sym;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Row-major SymmetricMatrix default constructor (HybridMatrix)";

      const blaze::SymmetricMatrix< blaze::HybridMatrix<VT,3UL,3UL,blaze::rowMajor> > sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Row-major SymmetricMatrix default constructor (DynamicMatrix)";

      const ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   // Size constructor (HybridMatrix)
   {
      test_ = "Row-major SymmetricMatrix size constructor (HybridMatrix)";

      const blaze::SymmetricMatrix< blaze::HybridMatrix<VT,3UL,3UL,blaze::rowMajor> > sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
   }

   // Size constructor (DynamicMatrix)
   {
      test_ = "Row-major SymmetricMatrix size constructor (DynamicMatrix)";

      const ST sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Row-major custom matrix constructors
   //=====================================================================================

   // Custom matrix constructor (ElementType*, size_t)
   {
      test_ = "Row-major SymmetricMatrix custom matrix constructor (ElementType*, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<VT,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<VT[]> array( new VT[5UL] );
      array[1] = vec( 1 );
      array[2] = vec( 2 );
      array[3] = vec( 2 );
      array[4] = vec( 1 );
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( array.get()+1UL, 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 1 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 1 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t)
   {
      test_ = "Row-major SymmetricMatrix custom matrix constructor (ElementType*, size_t, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<VT,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<VT[]> array( new VT[11UL] );
      array[1] = vec( 1 );
      array[2] = vec( 2 );
      array[6] = vec( 2 );
      array[7] = vec( 1 );
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( array.get()+1UL, 2UL, 5UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 1 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 1 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, Deleter)
   {
      test_ = "Row-major SymmetricMatrix custom matrix constructor (ElementType*, size_t, Deleter)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<VT,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<VT[]> array( new VT[4UL] );
      array[0] = vec( 1 );
      array[1] = vec( 2 );
      array[2] = vec( 2 );
      array[3] = vec( 1 );
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( array.release(), 2UL, blaze::ArrayDelete() );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 1 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 1 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t, Deleter)
   {
      test_ = "Row-major SymmetricMatrix custom matrix constructor (ElementType*, size_t, size_t, Deleter)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      typedef blaze::CustomMatrix<VT,unaligned,unpadded,rowMajor>  UnalignedUnpadded;
      std::unique_ptr<VT[]> array( new VT[10UL] );
      array[0] = vec( 1 );
      array[1] = vec( 2 );
      array[5] = vec( 2 );
      array[6] = vec( 1 );
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( array.release(), 2UL, 5UL, blaze::ArrayDelete() );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 1 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 1 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Row-major SymmetricMatrix copy constructor (0x0)";

      const ST sym1;
      const ST sym2( sym1 );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Row-major SymmetricMatrix copy constructor (3x3)";

      ST sym1( 3UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      const ST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Row-major SymmetricMatrix move constructor (0x0)";

      ST sym1;
      ST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Row-major SymmetricMatrix move constructor (3x3)";

      ST sym1( 3UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<VT,blaze::rowMajor> mat;
      const ST sym( mat );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Conversion constructor (symmetric)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      const ST sym( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-symmetric)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         const ST sym( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-symmetric SymmetricMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (SymmetricMatrix)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      const ST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Column-major SymmetricMatrix default constructor (StaticMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> > sym;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Column-major SymmetricMatrix default constructor (HybridMatrix)";

      blaze::SymmetricMatrix< blaze::HybridMatrix<VT,3UL,3UL,blaze::columnMajor> > sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Column-major SymmetricMatrix default constructor (DynamicMatrix)";

      OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   // Size constructor (HybridMatrix)
   {
      test_ = "Column-major SymmetricMatrix size constructor (HybridMatrix)";

      const blaze::SymmetricMatrix< blaze::HybridMatrix<VT,3UL,3UL,blaze::columnMajor> > sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkNonZeros( sym, 0UL );
   }

   // Size constructor (DynamicMatrix)
   {
      test_ = "Column-major SymmetricMatrix size constructor (DynamicMatrix)";

      const OST sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major custom matrix constructors
   //=====================================================================================

   // Custom matrix constructor (ElementType*, size_t)
   {
      test_ = "Column-major SymmetricMatrix custom matrix constructor (ElementType*, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<VT,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<VT[]> array( new VT[5UL] );
      array[1] = vec( 1 );
      array[2] = vec( 2 );
      array[3] = vec( 2 );
      array[4] = vec( 1 );
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( array.get()+1UL, 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 1 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 1 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t)
   {
      test_ = "Column-major SymmetricMatrix custom matrix constructor (ElementType*, size_t, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<VT,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<VT[]> array( new VT[11UL] );
      array[1] = vec( 1 );
      array[2] = vec( 2 );
      array[6] = vec( 2 );
      array[7] = vec( 1 );
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( array.get()+1UL, 2UL, 5UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 1 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 1 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, Deleter)
   {
      test_ = "Column-major SymmetricMatrix custom matrix constructor (ElementType*, size_t, Deleter)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<VT,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<VT[]> array( new VT[4UL] );
      array[0] = vec( 1 );
      array[1] = vec( 2 );
      array[2] = vec( 2 );
      array[3] = vec( 1 );
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( array.release(), 2UL, blaze::ArrayDelete() );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 1 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 1 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t, Deleter)
   {
      test_ = "Column-major SymmetricMatrix custom matrix constructor (ElementType*, size_t, size_t, Deleter)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      typedef blaze::CustomMatrix<VT,unaligned,unpadded,columnMajor>  UnalignedUnpadded;
      std::unique_ptr<VT[]> array( new VT[10UL] );
      array[0] = vec( 1 );
      array[1] = vec( 2 );
      array[5] = vec( 2 );
      array[6] = vec( 1 );
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( array.release(), 2UL, 5UL, blaze::ArrayDelete() );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 1 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 1 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Column-major SymmetricMatrix copy constructor (0x0)";

      const OST sym1;
      const OST sym2( sym1 );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Column-major SymmetricMatrix copy constructor (3x3)";

      OST sym1( 3UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      const OST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Column-major SymmetricMatrix move constructor (0x0)";

      OST sym1;
      OST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Column-major SymmetricMatrix move constructor (3x3)";

      OST sym1( 3UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<VT,blaze::columnMajor> mat;
      const OST sym( mat );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Conversion constructor (symmetric)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      const OST sym( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-symmetric)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         const OST sym( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-symmetric SymmetricMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (SymmetricMatrix)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      const OST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testAssignment()
{
   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix copy assignment (0x0)";

      ST sym1, sym2;

      sym2 = sym1;

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Row-major SymmetricMatrix copy assignment (3x3)";

      ST sym1( 3UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix move assignment (0x0)";

      ST sym1, sym2;

      sym2 = std::move( sym1 );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Row-major SymmetricMatrix move assignment (3x3)";

      ST sym1( 3UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = std::move( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<VT,blaze::rowMajor> mat;

      ST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Row-major/row-major dense matrix assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         ST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         ST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<VT,blaze::rowMajor> mat;

      ST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         ST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         ST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix computation assignment
   //=====================================================================================

   // Row-major/row-major dense matrix computation assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix computation assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix computation assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix computation assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix computation assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix computation assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         ST sym;
         sym = eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix computation assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix computation assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         ST sym;
         sym = eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix copy assignment (0x0)";

      OST sym1, sym2;

      sym2 = sym1;

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Column-major SymmetricMatrix copy assignment (3x3)";

      OST sym1( 3UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix move assignment (0x0)";

      OST sym1, sym2;

      sym2 = std::move( sym1 );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Column-major SymmetricMatrix move assignment (3x3)";

      OST sym1( 3UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = std::move( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<VT,blaze::columnMajor> mat;

      OST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Column-major/row-major dense matrix assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         OST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         OST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<VT,blaze::columnMajor> mat;

      OST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         OST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         OST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix computation assignment
   //=====================================================================================

   // Column-major/row-major dense matrix computation assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix computation assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix computation assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix computation assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix computation assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix computation assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         OST sym;
         sym = eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix computation assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix computation assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         OST sym;
         sym = eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
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
void DenseNonNumericTest::testAddAssign()
{
   //=====================================================================================
   // Row-major addition assignment
   //=====================================================================================

   // Row-major/row-major addition assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major addition assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major addition assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major addition assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major addition assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix addition assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major addition assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix addition assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major computation addition assignment
   //=====================================================================================

   // Row-major/row-major computation addition assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix computation addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major computation addition assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix computation addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major computation addition assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix computation addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major computation addition assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix computation addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major addition assignment
   //=====================================================================================

   // Column-major/row-major addition assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major addition assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major addition assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major addition assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major addition assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix addition assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major addition assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix addition assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major computation addition assignment
   //=====================================================================================

   // Column-major/row-major computation addition assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major computation addition assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major computation addition assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major computation addition assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
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
void DenseNonNumericTest::testSubAssign()
{
   //=====================================================================================
   // Row-major subtraction assignment
   //=====================================================================================

   // Row-major/row-major subtraction assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major subtraction assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major subtraction assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix subtraction assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec(  2 );
      sym1(0,2) = vec( -6 );
      sym1(1,1) = vec( -3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major subtraction assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix subtraction assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec(  2 );
      sym1(0,2) = vec( -6 );
      sym1(1,1) = vec( -3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major computation subtraction assignment
   //=====================================================================================

   // Row-major/row-major subtraction assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major subtraction assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major subtraction assignment
   //=====================================================================================

   // Column-major/row-major subtraction assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major subtraction assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n( ( -1 ) (  3 ) (   ) )\n( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major subtraction assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix subtraction assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec(  2 );
      sym1(0,2) = vec( -6 );
      sym1(1,1) = vec( -3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major subtraction assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix subtraction assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec(  2 );
      sym1(0,2) = vec( -6 );
      sym1(1,1) = vec( -3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major computation subtraction assignment
   //=====================================================================================

   // Column-major/row-major computation subtraction assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major computation subtraction assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major computation subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major computation subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
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
void DenseNonNumericTest::testMultAssign()
{
   //=====================================================================================
   // Row-major multiplication assignment
   //=====================================================================================

   // Row-major/row-major multiplication assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,0) = vec( 2 );
      mat(1,1) = vec( 2 );
      mat(2,2) = vec( 2 );

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec(  2 ) || sym(0,1) != vec( -8 ) || sym(0,2) != vec( 14 ) ||
          sym(1,0) != vec( -8 ) || sym(1,1) != vec(  4 ) || sym(1,2) != vec(  0 ) ||
          sym(2,0) != vec( 14 ) || sym(2,1) != vec(  0 ) || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major multiplication assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,0) = vec( 2 );
      mat(1,1) = vec( 2 );
      mat(2,2) = vec( 2 );

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec(  2 ) || sym(0,1) != vec( -8 ) || sym(0,2) != vec( 14 ) ||
          sym(1,0) != vec( -8 ) || sym(1,1) != vec(  4 ) || sym(1,2) != vec(  0 ) ||
          sym(2,0) != vec( 14 ) || sym(2,1) != vec(  0 ) || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major multiplication assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major multiplication assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major multiplication assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix multiplication assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,0) = vec( 2 );
      sym1(0,1) = vec( 0 );
      sym1(0,2) = vec( 0 );
      sym1(1,1) = vec( 2 );
      sym1(1,2) = vec( 0 );
      sym1(2,2) = vec( 2 );

      ST sym2( 3UL );
      sym2(0,0) = vec(  1 );
      sym2(0,1) = vec( -4 );
      sym2(0,2) = vec(  7 );
      sym2(1,1) = vec(  2 );
      sym2(1,2) = vec(  0 );
      sym2(2,2) = vec(  3 );

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 9UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 3UL );
      checkNonZeros( sym2, 2UL, 3UL );

      if( sym2(0,0) != vec(  2 ) || sym2(0,1) != vec( -8 ) || sym2(0,2) != vec( 14 ) ||
          sym2(1,0) != vec( -8 ) || sym2(1,1) != vec(  4 ) || sym2(1,2) != vec(  0 ) ||
          sym2(2,0) != vec( 14 ) || sym2(2,1) != vec(  0 ) || sym2(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major multiplication assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix multiplication assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,0) = vec( 2 );
      sym1(0,1) = vec( 0 );
      sym1(0,2) = vec( 0 );
      sym1(1,1) = vec( 2 );
      sym1(1,2) = vec( 0 );
      sym1(2,2) = vec( 2 );

      ST sym2( 3UL );
      sym2(0,0) = vec(  1 );
      sym2(0,1) = vec( -4 );
      sym2(0,2) = vec(  7 );
      sym2(1,1) = vec(  2 );
      sym2(1,2) = vec(  0 );
      sym2(2,2) = vec(  3 );

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 9UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 3UL );
      checkNonZeros( sym2, 2UL, 3UL );

      if( sym2(0,0) != vec(  2 ) || sym2(0,1) != vec( -8 ) || sym2(0,2) != vec( 14 ) ||
          sym2(1,0) != vec( -8 ) || sym2(1,1) != vec(  4 ) || sym2(1,2) != vec(  0 ) ||
          sym2(2,0) != vec( 14 ) || sym2(2,1) != vec(  0 ) || sym2(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major multiplication assignment
   //=====================================================================================

   // Column-major/row-major multiplication assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,0) = vec( 2 );
      mat(1,1) = vec( 2 );
      mat(2,2) = vec( 2 );

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec(  2 ) || sym(0,1) != vec( -8 ) || sym(0,2) != vec( 14 ) ||
          sym(1,0) != vec( -8 ) || sym(1,1) != vec(  4 ) || sym(1,2) != vec(  0 ) ||
          sym(2,0) != vec( 14 ) || sym(2,1) != vec(  0 ) || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major multiplication assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,0) = vec( 2 );
      mat(1,1) = vec( 2 );
      mat(2,2) = vec( 2 );

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec(  2 ) || sym(0,1) != vec( -8 ) || sym(0,2) != vec( 14 ) ||
          sym(1,0) != vec( -8 ) || sym(1,1) != vec(  4 ) || sym(1,2) != vec(  0 ) ||
          sym(2,0) != vec( 14 ) || sym(2,1) != vec(  0 ) || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major multiplication assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major multiplication assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major multiplication assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix multiplication assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,0) = vec( 2 );
      sym1(0,1) = vec( 0 );
      sym1(0,2) = vec( 0 );
      sym1(1,1) = vec( 2 );
      sym1(1,2) = vec( 0 );
      sym1(2,2) = vec( 2 );

      OST sym2( 3UL );
      sym2(0,0) = vec(  1 );
      sym2(0,1) = vec( -4 );
      sym2(0,2) = vec(  7 );
      sym2(1,1) = vec(  2 );
      sym2(1,2) = vec(  0 );
      sym2(2,2) = vec(  3 );

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 9UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 3UL );
      checkNonZeros( sym2, 2UL, 3UL );

      if( sym2(0,0) != vec(  2 ) || sym2(0,1) != vec( -8 ) || sym2(0,2) != vec( 14 ) ||
          sym2(1,0) != vec( -8 ) || sym2(1,1) != vec(  4 ) || sym2(1,2) != vec(  0 ) ||
          sym2(2,0) != vec( 14 ) || sym2(2,1) != vec(  0 ) || sym2(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major multiplication assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix multiplication assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,0) = vec( 2 );
      sym1(0,1) = vec( 0 );
      sym1(0,2) = vec( 0 );
      sym1(1,1) = vec( 2 );
      sym1(1,2) = vec( 0 );
      sym1(2,2) = vec( 2 );

      OST sym2( 3UL );
      sym2(0,0) = vec(  1 );
      sym2(0,1) = vec( -4 );
      sym2(0,2) = vec(  7 );
      sym2(1,1) = vec(  2 );
      sym2(1,2) = vec(  0 );
      sym2(2,2) = vec(  3 );

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 9UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 3UL );
      checkNonZeros( sym2, 2UL, 3UL );

      if( sym2(0,0) != vec(  2 ) || sym2(0,1) != vec( -8 ) || sym2(0,2) != vec( 14 ) ||
          sym2(1,0) != vec( -8 ) || sym2(1,1) != vec(  4 ) || sym2(1,2) != vec(  0 ) ||
          sym2(2,0) != vec( 14 ) || sym2(2,1) != vec(  0 ) || sym2(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      ST sym( 3UL );
      sym(1,2) = vec(  1 );
      sym(2,0) = vec( -2 );
      sym(2,2) = vec(  3 );

      sym *= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -4 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  2 ) ||
          sym(2,0) != vec( -4 )  || sym(2,1) != vec( 2 )   || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -4 ) )\n"
                                     "( (    ) (   ) (  2 ) )\n"
                                     "( ( -4 ) ( 2 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      ST sym( 3UL );
      sym(1,2) = vec(  1 );
      sym(2,0) = vec( -2 );
      sym(2,2) = vec(  3 );

      sym = sym * 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -4 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  2 ) ||
          sym(2,0) != vec( -4 )  || sym(2,1) != vec( 2 )   || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -4 ) )\n"
                                     "( (    ) (   ) (  2 ) )\n"
                                     "( ( -4 ) ( 2 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      ST sym( 3UL );
      sym(1,2) = vec(  1 );
      sym(2,0) = vec( -2 );
      sym(2,2) = vec(  3 );

      sym = 2 * sym;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -4 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  2 ) ||
          sym(2,0) != vec( -4 )  || sym(2,1) != vec( 2 )   || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -4 ) )\n"
                                     "( (    ) (   ) (  2 ) )\n"
                                     "( ( -4 ) ( 2 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      ST sym( 3UL );
      sym(1,2) = vec(  2 );
      sym(2,0) = vec( -4 );
      sym(2,2) = vec(  6 );

      sym /= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -2 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  1 ) ||
          sym(2,0) != vec( -2 )  || sym(2,1) != vec( 1 )   || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -2 ) )\n"
                                     "( (    ) (   ) (  1 ) )\n"
                                     "( ( -2 ) ( 1 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      ST sym( 3UL );
      sym(1,2) = vec(  2 );
      sym(2,0) = vec( -4 );
      sym(2,2) = vec(  6 );

      sym = sym / 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -2 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  1 ) ||
          sym(2,0) != vec( -2 )  || sym(2,1) != vec( 1 )   || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -2 ) )\n"
                                     "( (    ) (   ) (  1 ) )\n"
                                     "( ( -2 ) ( 1 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major SymmetricMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::scale()";

      // Initialization check
      ST sym( 3UL );
      sym(1,2) = vec(  1 );
      sym(2,0) = vec( -2 );
      sym(2,2) = vec(  3 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -2 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  1 ) ||
          sym(2,0) != vec( -2 )  || sym(2,1) != vec( 1 )   || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -2 ) )\n"
                                     "( (    ) (   ) (  1 ) )\n"
                                     "( ( -2 ) ( 1 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      sym.scale( 2 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -4 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  2 ) ||
          sym(2,0) != vec( -4 )  || sym(2,1) != vec( 2 )   || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -4 ) )\n"
                                     "( (    ) (   ) (  2 ) )\n"
                                     "( ( -4 ) ( 2 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      sym.scale( 0.5 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -2 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  1 ) ||
          sym(2,0) != vec( -2 )  || sym(2,1) != vec( 1 )   || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -2 ) )\n"
                                     "( (    ) (   ) (  1 ) )\n"
                                     "( ( -2 ) ( 1 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      OST sym( 3UL );
      sym(1,2) = vec(  1 );
      sym(2,0) = vec( -2 );
      sym(2,2) = vec(  3 );

      sym *= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -4 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  2 ) ||
          sym(2,0) != vec( -4 )  || sym(2,1) != vec( 2 )   || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -4 ) )\n"
                                     "( (    ) (   ) (  2 ) )\n"
                                     "( ( -4 ) ( 2 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      OST sym( 3UL );
      sym(1,2) = vec(  1 );
      sym(2,0) = vec( -2 );
      sym(2,2) = vec(  3 );

      sym = sym * 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -4 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  2 ) ||
          sym(2,0) != vec( -4 )  || sym(2,1) != vec( 2 )   || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -4 ) )\n"
                                     "( (    ) (   ) (  2 ) )\n"
                                     "( ( -4 ) ( 2 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      OST sym( 3UL );
      sym(1,2) = vec(  1 );
      sym(2,0) = vec( -2 );
      sym(2,2) = vec(  3 );

      sym = 2 * sym;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -4 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  2 ) ||
          sym(2,0) != vec( -4 )  || sym(2,1) != vec( 2 )   || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -4 ) )\n"
                                     "( (    ) (   ) (  2 ) )\n"
                                     "( ( -4 ) ( 2 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      OST sym( 3UL );
      sym(1,2) = vec(  2 );
      sym(2,0) = vec( -4 );
      sym(2,2) = vec(  6 );

      sym /= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -2 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  1 ) ||
          sym(2,0) != vec( -2 )  || sym(2,1) != vec( 1 )   || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -2 ) )\n"
                                     "( (    ) (   ) (  1 ) )\n"
                                     "( ( -2 ) ( 1 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      OST sym( 3UL );
      sym(1,2) = vec(  2 );
      sym(2,0) = vec( -4 );
      sym(2,2) = vec(  6 );

      sym = sym / 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -2 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  1 ) ||
          sym(2,0) != vec( -2 )  || sym(2,1) != vec( 1 )   || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -2 ) )\n"
                                     "( (    ) (   ) (  1 ) )\n"
                                     "( ( -2 ) ( 1 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major SymmetricMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::scale()";

      // Initialization check
      OST sym( 3UL );
      sym(1,2) = vec(  1 );
      sym(2,0) = vec( -2 );
      sym(2,2) = vec(  3 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -2 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  1 ) ||
          sym(2,0) != vec( -2 )  || sym(2,1) != vec( 1 )   || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -2 ) )\n"
                                     "( (    ) (   ) (  1 ) )\n"
                                     "( ( -2 ) ( 1 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      sym.scale( 2 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -4 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  2 ) ||
          sym(2,0) != vec( -4 )  || sym(2,1) != vec( 2 )   || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -4 ) )\n"
                                     "( (    ) (   ) (  2 ) )\n"
                                     "( ( -4 ) ( 2 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      sym.scale( 0.5 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( -2 ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || sym(1,2) != vec(  1 ) ||
          sym(2,0) != vec( -2 )  || sym(2,1) != vec( 1 )   || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) (   ) ( -2 ) )\n"
                                     "( (    ) (   ) (  1 ) )\n"
                                     "( ( -2 ) ( 1 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::operator()";

      ST sym( 3UL );

      // Writing the element (1,1)
      sym(1,1) = vec( 1 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 1UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 0UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || !isDefault( sym(0,2) ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || !isDefault( sym(1,2) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) (   ) )\n"
                                     "( (   ) ( 1 ) (   ) )\n"
                                     "( (   ) (   ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the elements (2,1) and (1,2)
      sym(2,1) = vec( 2 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || !isDefault( sym(0,2) ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 2 )   ||
          !isDefault( sym(2,0) ) || sym(2,1) != vec( 2 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) (   ) )\n"
                                     "( (   ) ( 1 ) ( 2 ) )\n"
                                     "( (   ) ( 2 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the elements (0,2) and (2,0)
      sym(0,2) = sym(1,2);

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 2 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 2 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 2 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Adding to the elements (1,2) and (2,1)
      sym(1,2) += vec( 3 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 5 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 1 ) ( 5 ) )\n"
                                     "( ( 2 ) ( 5 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtracting from the elements (0,1) and (1,0)
      sym(1,2) -= vec( 4 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 1 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 1 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 1 ) ( 1 ) )\n"
                                     "( ( 2 ) ( 1 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplying the element (1,1)
      sym(2,0) *= 3;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 6 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 1 ) ||
          sym(2,0) != vec( 6 )   || sym(2,1) != vec( 1 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 6 ) )\n"
                                     "( (   ) ( 1 ) ( 1 ) )\n"
                                     "( ( 6 ) ( 1 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Dividing the elements (0,2) and (2,0)
      sym(2,0) /= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 3 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 1 ) ||
          sym(2,0) != vec( 3 )   || sym(2,1) != vec( 1 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 3 ) )\n"
                                     "( (   ) ( 1 ) ( 1 ) )\n"
                                     "( ( 3 ) ( 1 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::operator()";

      OST sym( 3UL );

      // Writing the element (1,1)
      sym(1,1) = vec( 1 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 1UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 0UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || !isDefault( sym(0,2) ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || !isDefault( sym(1,2) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) (   ) )\n"
                                     "( (   ) ( 1 ) (   ) )\n"
                                     "( (   ) (   ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the elements (2,1) and (1,2)
      sym(2,1) = vec( 2 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || !isDefault( sym(0,2) ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 2 )   ||
          !isDefault( sym(2,0) ) || sym(2,1) != vec( 2 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) (   ) )\n"
                                     "( (   ) ( 1 ) ( 2 ) )\n"
                                     "( (   ) ( 2 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the elements (0,2) and (2,0)
      sym(0,2) = sym(1,2);

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 2 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 2 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 2 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Adding to the elements (1,2) and (2,1)
      sym(1,2) += vec( 3 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 5 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 1 ) ( 5 ) )\n"
                                     "( ( 2 ) ( 5 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtracting from the elements (0,1) and (1,0)
      sym(1,2) -= vec( 4 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 1 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 1 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 1 ) ( 1 ) )\n"
                                     "( ( 2 ) ( 1 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplying the element (1,1)
      sym(2,0) *= 3;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 6 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 1 ) ||
          sym(2,0) != vec( 6 )   || sym(2,1) != vec( 1 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 6 ) )\n"
                                     "( (   ) ( 1 ) ( 1 ) )\n"
                                     "( ( 6 ) ( 1 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Dividing the elements (0,2) and (2,0)
      sym(2,0) /= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || sym(0,2) != vec( 3 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 1 )   || sym(1,2) != vec( 1 ) ||
          sym(2,0) != vec( 3 )   || sym(2,1) != vec( 1 )   || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) ( 3 ) )\n"
                                     "( (   ) ( 1 ) ( 1 ) )\n"
                                     "( ( 3 ) ( 1 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      typedef ST::Iterator       Iterator;
      typedef ST::ConstIterator  ConstIterator;

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(1,2) = vec( 2 );
      sym(2,2) = vec( 3 );

      // Testing the Iterator default constructor
      {
         test_ = "Row-major Iterator default constructor";

         Iterator it = Iterator();

         if( it != Iterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing the ConstIterator default constructor
      {
         test_ = "Row-major ConstIterator default constructor";

         ConstIterator it = ConstIterator();

         if( it != ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Row-major Iterator/ConstIterator conversion";

         ConstIterator it( begin( sym, 1UL ) );

         if( it == end( sym, 1UL ) || *it != vec( 1 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator
      {
         test_ = "Row-major Iterator subtraction";

         const size_t number( end( sym, 0UL ) - begin( sym, 0UL ) );

         if( number != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st row via ConstIterator
      {
         test_ = "Row-major ConstIterator subtraction";

         const size_t number( cend( sym, 1UL ) - cbegin( sym, 1UL ) );

         if( number != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing read-only access via ConstIterator
      {
         test_ = "Row-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( sym, 2UL ) );
         ConstIterator end( cend( sym, 2UL ) );

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != vec( 2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != vec( 2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != vec( 3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 2UL;

         if( it == end || *it != vec( 3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 2UL;

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar subtraction failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = 3UL + it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Scalar/iterator addition failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Row-major assignment via Iterator";

         int value = 7;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it = vec( value++ );
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 )   || sym(0,2) != vec( 7 ) ||
             sym(1,0) != vec( 1 )   || !isDefault( sym(1,1) ) || sym(1,2) != vec( 8 ) ||
             sym(2,0) != vec( 7 )   || sym(2,1) != vec( 8 )   || sym(2,2) != vec( 9 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (   ) ( 1 ) ( 7 ) )\n"
                                        "( ( 1 ) (   ) ( 8 ) )\n"
                                        "( ( 7 ) ( 8 ) ( 9 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it += vec( value++ );
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec(  1 )  || sym(0,2) != vec( 11 ) ||
             sym(1,0) != vec(  1 )  || !isDefault( sym(1,1) ) || sym(1,2) != vec( 13 ) ||
             sym(2,0) != vec( 11 )  || sym(2,1) != vec( 13 )  || sym(2,2) != vec( 15 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (    ) (  1 ) ( 11 ) )\n"
                                        "( (  1 ) (    ) ( 13 ) )\n"
                                        "( ( 11 ) ( 13 ) ( 15 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it -= vec( value++ );
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 )   || sym(0,2) != vec( 7 ) ||
             sym(1,0) != vec( 1 )   || !isDefault( sym(1,1) ) || sym(1,2) != vec( 8 ) ||
             sym(2,0) != vec( 7 )   || sym(2,1) != vec( 8 )   || sym(2,2) != vec( 9 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (   ) ( 1 ) ( 7 ) )\n"
                                        "( ( 1 ) (   ) ( 8 ) )\n"
                                        "( ( 7 ) ( 8 ) ( 9 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         int value = 2;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it *= value++;
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec(  1 )  || sym(0,2) != vec( 14 ) ||
             sym(1,0) != vec(  1 )  || !isDefault( sym(1,1) ) || sym(1,2) != vec( 24 ) ||
             sym(2,0) != vec( 14 )  || sym(2,1) != vec( 24 )  || sym(2,2) != vec( 36 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (    ) (  1 ) ( 14 ) )\n"
                                        "( (  1 ) (    ) ( 24 ) )\n"
                                        "( ( 14 ) ( 24 ) ( 36 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it /= 2;
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec(  1 )   || sym(0,2) != vec(  7 ) ||
             sym(1,0) != vec( 1 )   || !isDefault( sym(1,1)  ) || sym(1,2) != vec( 12 ) ||
             sym(2,0) != vec( 7 )   || sym(2,1) != vec( 12 )   || sym(2,2) != vec( 18 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (   ) (  1 ) (  7 ) )\n"
                                        "( ( 1 ) (    ) ( 12 ) )\n"
                                        "( ( 7 ) ( 12 ) ( 18 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      typedef OST::Iterator       Iterator;
      typedef OST::ConstIterator  ConstIterator;

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(1,2) = vec( 2 );
      sym(2,2) = vec( 3 );

      // Testing the Iterator default constructor
      {
         test_ = "Column-major Iterator default constructor";

         Iterator it = Iterator();

         if( it != Iterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing the ConstIterator default constructor
      {
         test_ = "Column-major ConstIterator default constructor";

         ConstIterator it = ConstIterator();

         if( it != ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Column-major Iterator/ConstIterator conversion";

         ConstIterator it( begin( sym, 1UL ) );

         if( it == end( sym, 1UL ) || *it != vec( 1 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator
      {
         test_ = "Column-major Iterator subtraction";

         const size_t number( end( sym, 0UL ) - begin( sym, 0UL ) );

         if( number != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st column via ConstIterator
      {
         test_ = "Column-major ConstIterator subtraction";

         const size_t number( cend( sym, 1UL ) - cbegin( sym, 1UL ) );

         if( number != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing read-only access via ConstIterator
      {
         test_ = "Column-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( sym, 2UL ) );
         ConstIterator end( cend( sym, 2UL ) );

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != vec( 2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != vec( 2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != vec( 3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 2UL;

         if( it == end || *it != vec( 3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 2UL;

         if( it == end || !isDefault( *it ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar subtraction failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = 3UL + it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Scalar/iterator addition failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Column-major assignment via Iterator";

         int value = 7;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it = vec( value++ );
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 )   || sym(0,2) != vec( 7 ) ||
             sym(1,0) != vec( 1 )   || !isDefault( sym(1,1) ) || sym(1,2) != vec( 8 ) ||
             sym(2,0) != vec( 7 )   || sym(2,1) != vec( 8 )   || sym(2,2) != vec( 9 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (   ) ( 1 ) ( 7 ) )\n"
                                        "( ( 1 ) (   ) ( 8 ) )\n"
                                        "( ( 7 ) ( 8 ) ( 9 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it += vec( value++ );
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec(  1 )  || sym(0,2) != vec( 11 ) ||
             sym(1,0) != vec(  1 )  || !isDefault( sym(1,1) ) || sym(1,2) != vec( 13 ) ||
             sym(2,0) != vec( 11 )  || sym(2,1) != vec( 13 )  || sym(2,2) != vec( 15 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (    ) (  1 ) ( 11 ) )\n"
                                        "( (  1 ) (    ) ( 13 ) )\n"
                                        "( ( 11 ) ( 13 ) ( 15 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it -= vec( value++ );
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 )   || sym(0,2) != vec( 7 ) ||
             sym(1,0) != vec( 1 )   || !isDefault( sym(1,1) ) || sym(1,2) != vec( 8 ) ||
             sym(2,0) != vec( 7 )   || sym(2,1) != vec( 8 )   || sym(2,2) != vec( 9 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (   ) ( 1 ) ( 7 ) )\n"
                                        "( ( 1 ) (   ) ( 8 ) )\n"
                                        "( ( 7 ) ( 8 ) ( 9 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         int value = 2;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it *= value++;
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec(  1 )  || sym(0,2) != vec( 14 ) ||
             sym(1,0) != vec(  1 )  || !isDefault( sym(1,1) ) || sym(1,2) != vec( 24 ) ||
             sym(2,0) != vec( 14 )  || sym(2,1) != vec( 24 )  || sym(2,2) != vec( 36 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (    ) (  1 ) ( 14 ) )\n"
                                        "( (  1 ) (    ) ( 24 ) )\n"
                                        "( ( 14 ) ( 24 ) ( 36 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it /= 2;
         }

         if( !isDefault( sym(0,0) ) || sym(0,1) != vec(  1 )  || sym(0,2) != vec(  7 ) ||
             sym(1,0) != vec( 1 )   || !isDefault( sym(1,1) ) || sym(1,2) != vec( 12 ) ||
             sym(2,0) != vec( 7 )   || sym(2,1) != vec( 12 )  || sym(2,2) != vec( 18 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (   ) (  1 ) (  7 ) )\n"
                                        "( ( 1 ) (    ) ( 12 ) )\n"
                                        "( ( 7 ) ( 12 ) ( 18 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
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
void DenseNonNumericTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::nonZeros()";

      // Empty matrix
      {
         ST sym( 3UL );

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 0UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 0UL );
         checkNonZeros( sym, 2UL, 0UL );
      }

      // Partially filled matrix
      {
         ST sym( 3UL );
         sym(0,0) = vec(  2 );
         sym(1,2) = vec(  4 );
         sym(2,0) = VT();
         sym(2,2) = vec( -6 );

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 4UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
      }

      // Fully filled matrix
      {
         ST sym( 3UL );
         sym(0,0) = vec(  2 );
         sym(0,1) = vec( -4 );
         sym(0,2) = vec( -6 );
         sym(1,1) = vec(  8 );
         sym(1,2) = vec( 10 );
         sym(2,2) = vec( 12 );

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 3UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 3UL );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::nonZeros()";

      // Empty matrix
      {
         OST sym( 3UL );

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 0UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 0UL );
         checkNonZeros( sym, 2UL, 0UL );
      }

      // Partially filled matrix
      {
         OST sym( 3UL );
         sym(0,0) = vec(  2 );
         sym(1,2) = vec(  4 );
         sym(2,0) = VT();
         sym(2,2) = vec( -6 );

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 4UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
      }

      // Fully filled matrix
      {
         OST sym( 3UL );
         sym(0,0) = vec(  2 );
         sym(0,1) = vec( -4 );
         sym(0,2) = vec( -6 );
         sym(1,1) = vec(  8 );
         sym(1,2) = vec( 10 );
         sym(2,2) = vec( 12 );

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 3UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 3UL );
      }
   }
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
void DenseNonNumericTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::reset()";

      // Initialization check
      ST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,1) = vec( 2 );
      sym(0,2) = vec( 3 );
      sym(1,1) = vec( 4 );
      sym(1,2) = vec( 5 );
      sym(2,2) = vec( 6 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) || sym(0,2) != vec( 3 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 4 ) || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 3 ) || sym(2,1) != vec( 5 ) || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) ( 3 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n"
                                     "( ( 3 ) ( 5 ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( sym(0,1) );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 0 ) || sym(0,2) != vec( 3 ) ||
          sym(1,0) != vec( 0 ) || sym(1,1) != vec( 4 ) || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 3 ) || sym(2,1) != vec( 5 ) || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 0 ) ( 3 ) )\n"
                                     "( ( 0 ) ( 4 ) ( 5 ) )\n"
                                     "( ( 3 ) ( 5 ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( sym, 1UL );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 3 )   ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) )   || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 3 )   || !isDefault( sym(2,1) )   || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 3 ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( ( 3 ) (   ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 0UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || !isDefault( sym(0,2) ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::reset()";

      // Initialization check
      OST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,1) = vec( 2 );
      sym(0,2) = vec( 3 );
      sym(1,1) = vec( 4 );
      sym(1,2) = vec( 5 );
      sym(2,2) = vec( 6 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) || sym(0,2) != vec( 3 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 4 ) || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 3 ) || sym(2,1) != vec( 5 ) || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) ( 3 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n"
                                     "( ( 3 ) ( 5 ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( sym(0,1) );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 0 ) || sym(0,2) != vec( 3 ) ||
          sym(1,0) != vec( 0 ) || sym(1,1) != vec( 4 ) || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 3 ) || sym(2,1) != vec( 5 ) || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 0 ) ( 3 ) )\n"
                                     "( ( 0 ) ( 4 ) ( 5 ) )\n"
                                     "( ( 3 ) ( 5 ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( sym, 1UL );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 3 )   ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 3 )   || !isDefault( sym(2,1) ) || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 3 ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( ( 3 ) (   ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 0UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) || !isDefault( sym(0,2) ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::clear()";

      // Initialization check
      ST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,1) = vec( 2 );
      sym(0,2) = vec( 3 );
      sym(1,1) = vec( 4 );
      sym(1,2) = vec( 5 );
      sym(2,2) = vec( 6 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) || sym(0,2) != vec( 3 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 4 ) || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 3 ) || sym(2,1) != vec( 5 ) || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) ( 3 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n"
                                     "( ( 3 ) ( 5 ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( sym(0,1) );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 3 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 4 )   || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 3 )   || sym(2,1) != vec( 5 )   || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 3 ) )\n"
                                     "( (   ) ( 4 ) ( 5 ) )\n"
                                     "( ( 3 ) ( 5 ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( sym );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::clear()";

      // Initialization check
      OST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,1) = vec( 2 );
      sym(0,2) = vec( 3 );
      sym(1,1) = vec( 4 );
      sym(1,2) = vec( 5 );
      sym(2,2) = vec( 6 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) || sym(0,2) != vec( 3 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 4 ) || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 3 ) || sym(2,1) != vec( 5 ) || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) ( 3 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n"
                                     "( ( 3 ) ( 5 ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( sym(0,1) );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 3 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 4 )   || sym(1,2) != vec( 5 ) ||
          sym(2,0) != vec( 3 )   || sym(2,1) != vec( 5 )   || sym(2,2) != vec( 6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 3 ) )\n"
                                     "( (   ) ( 4 ) ( 5 ) )\n"
                                     "( ( 3 ) ( 5 ) ( 6 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( sym );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }
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
void DenseNonNumericTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::resize()";

      // Initialization check
      ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Resizing to 2x2
      sym.resize( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( ) ( ) )\n"
                                     "( ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      sym(0,1) = vec( 1 );
      sym(1,1) = vec( 2 );
      sym.resize( 4UL, true );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 16UL );
      checkNonZeros( sym,  3UL );
      checkNonZeros( sym,  0UL, 1UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 0UL );
      checkNonZeros( sym,  3UL, 0UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 )   || !isDefault( sym(0,2) ) || !isDefault( sym(0,3) ) ||
          sym(1,0) != vec( 1 )   || sym(1,1) != vec( 2 )   || !isDefault( sym(1,2) ) || !isDefault( sym(1,3) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) || !isDefault( sym(2,3) ) ||
          !isDefault( sym(3,0) ) || !isDefault( sym(3,1) ) || !isDefault( sym(3,2) ) || !isDefault( sym(3,3) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) ( 1 ) (   ) (   ) )\n"
                                     "( ( 1 ) ( 2 ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      sym(2,2) = vec( 3 );
      sym.resize( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 ) ||
          sym(1,0) != vec( 1 )   || sym(1,1) != vec( 2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) ( 1 ) )\n"
                                     "( ( 1 ) ( 2 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      sym.resize( 0UL );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::resize()";

      // Initialization check
      OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Resizing to 2x2
      sym.resize( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( ) ( ) )\n"
                                     "( ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      sym(0,1) = vec( 1 );
      sym(1,1) = vec( 2 );
      sym.resize( 4UL, true );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 16UL );
      checkNonZeros( sym,  3UL );
      checkNonZeros( sym,  0UL, 1UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 0UL );
      checkNonZeros( sym,  3UL, 0UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 )   || !isDefault( sym(0,2) ) || !isDefault( sym(0,3) ) ||
          sym(1,0) != vec( 1 )   || sym(1,1) != vec( 2 )   || !isDefault( sym(1,2) ) || !isDefault( sym(1,3) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) || !isDefault( sym(2,3) ) ||
          !isDefault( sym(3,0) ) || !isDefault( sym(3,1) ) || !isDefault( sym(3,2) ) || !isDefault( sym(3,3) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) ( 1 ) (   ) (   ) )\n"
                                     "( ( 1 ) ( 2 ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      sym(2,2) = vec( 3 );
      sym.resize( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 ) ||
          sym(1,0) != vec( 1 )   || sym(1,1) != vec( 2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) ( 1 ) )\n"
                                     "( ( 1 ) ( 2 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      sym.resize( 0UL );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c extend() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseNonNumericTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::extend()";

      // Initialization check
      ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Extending the size of the matrix to 2x2
      sym.extend( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( ) ( ) )\n"
                                     "( ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      sym(0,1) = vec( 1 );
      sym(1,1) = vec( 2 );
      sym.extend( 2UL, true );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 16UL );
      checkNonZeros( sym,  3UL );
      checkNonZeros( sym,  0UL, 1UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 0UL );
      checkNonZeros( sym,  3UL, 0UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 )   || !isDefault( sym(0,2) ) || !isDefault( sym(0,3) ) ||
          sym(1,0) != vec( 1 )   || sym(1,1) != vec( 2 )   || !isDefault( sym(1,2) ) || !isDefault( sym(1,3) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) || !isDefault( sym(2,3) ) ||
          !isDefault( sym(3,0) ) || !isDefault( sym(3,1) ) || !isDefault( sym(3,2) ) || !isDefault( sym(3,3) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) ( 1 ) (   ) (   ) )\n"
                                     "( ( 1 ) ( 2 ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::extend()";

      // Initialization check
      OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Extending the size of the matrix to 2x2
      sym.extend( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );

      if( !isDefault( sym(0,0) ) || !isDefault( sym(0,1) ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( ) ( ) )\n"
                                     "( ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      sym(0,1) = vec( 1 );
      sym(1,1) = vec( 2 );
      sym.extend( 2UL, true );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 16UL );
      checkNonZeros( sym,  3UL );
      checkNonZeros( sym,  0UL, 1UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 0UL );
      checkNonZeros( sym,  3UL, 0UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( 1 )   || !isDefault( sym(0,2) ) || !isDefault( sym(0,3) ) ||
          sym(1,0) != vec( 1 )   || sym(1,1) != vec( 2 )   || !isDefault( sym(1,2) ) || !isDefault( sym(1,3) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) || !isDefault( sym(2,3) ) ||
          !isDefault( sym(3,0) ) || !isDefault( sym(3,1) ) || !isDefault( sym(3,2) ) || !isDefault( sym(3,3) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (   ) ( 1 ) (   ) (   ) )\n"
                                     "( ( 1 ) ( 2 ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::reserve()";

      // Initialization check
      ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Increasing the capacity of the matrix
      sym.reserve( 10UL );

      checkRows    ( sym,  0UL );
      checkColumns ( sym,  0UL );
      checkCapacity( sym, 10UL );
      checkNonZeros( sym,  0UL );

      // Further increasing the capacity of the matrix
      sym.reserve( 20UL );

      checkRows    ( sym,  0UL );
      checkColumns ( sym,  0UL );
      checkCapacity( sym, 20UL );
      checkNonZeros( sym,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::reserve()";

      // Initialization check
      OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Increasing the capacity of the matrix
      sym.reserve( 10UL );

      checkRows    ( sym,  0UL );
      checkColumns ( sym,  0UL );
      checkCapacity( sym, 10UL );
      checkNonZeros( sym,  0UL );

      // Further increasing the capacity of the matrix
      sym.reserve( 20UL );

      checkRows    ( sym,  0UL );
      checkColumns ( sym,  0UL );
      checkCapacity( sym, 20UL );
      checkNonZeros( sym,  0UL );
   }
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
void DenseNonNumericTest::testTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via transpose()";

      ST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 3 );
      sym(1,2) = vec( 4 );
      sym(2,2) = vec( 5 );

      transpose( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 3 )   || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 4 )   || sym(2,2) != vec( 5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 3 ) ( 4 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via trans()";

      ST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 3 );
      sym(1,2) = vec( 4 );
      sym(2,2) = vec( 5 );

      sym = trans( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 3 )   || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 4 )   || sym(2,2) != vec( 5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 3 ) ( 4 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      OST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 3 );
      sym(1,2) = vec( 4 );
      sym(2,2) = vec( 5 );

      transpose( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 3 )   || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 4 )   || sym(2,2) != vec( 5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 3 ) ( 4 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via trans()";

      OST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 3 );
      sym(1,2) = vec( 4 );
      sym(2,2) = vec( 5 );

      sym = trans( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 3 )   || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 4 )   || sym(2,2) != vec( 5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 3 ) ( 4 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the SymmetricMatrix
// specialization. Additionally, it performs a test of self-transpose via the \c ctrans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseNonNumericTest::testCTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via ctranspose()";

      ST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 3 );
      sym(1,2) = vec( 4 );
      sym(2,2) = vec( 5 );

      ctranspose( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 3 )   || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 4 )   || sym(2,2) != vec( 5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 3 ) ( 4 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via ctrans()";

      ST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 3 );
      sym(1,2) = vec( 4 );
      sym(2,2) = vec( 5 );

      sym = ctrans( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 3 )   || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 4 )   || sym(2,2) != vec( 5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 3 ) ( 4 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      OST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 3 );
      sym(1,2) = vec( 4 );
      sym(2,2) = vec( 5 );

      ctranspose( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 3 )   || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 4 )   || sym(2,2) != vec( 5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 3 ) ( 4 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via ctrans()";

      OST sym( 3UL );
      sym(0,0) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 3 );
      sym(1,2) = vec( 4 );
      sym(2,2) = vec( 5 );

      sym = ctrans( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 2 ) ||
          !isDefault( sym(1,0) ) || sym(1,1) != vec( 3 )   || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 2 )   || sym(2,1) != vec( 4 )   || sym(2,2) != vec( 5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 2 ) )\n"
                                     "( (   ) ( 3 ) ( 4 ) )\n"
                                     "( ( 2 ) ( 4 ) ( 5 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix swap";

      ST sym1( 2UL );
      sym1(0,0) = vec( 1 );
      sym1(0,1) = vec( 2 );
      sym1(1,1) = vec( 3 );

      ST sym2( 2UL );
      sym2(0,0) = vec( 4 );
      sym2(0,1) = vec( 5 );

      swap( sym1, sym2 );

      checkRows    ( sym1, 2UL );
      checkColumns ( sym1, 2UL );
      checkCapacity( sym1, 4UL );
      checkNonZeros( sym1, 3UL );
      checkNonZeros( sym1, 0UL, 2UL );
      checkNonZeros( sym1, 1UL, 1UL );

      if( sym1(0,0) != vec( 4 ) || sym1(0,1) != vec( 5 ) ||
          sym1(1,0) != vec( 5 ) || !isDefault( sym1(1,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym1 << "\n"
             << "   Expected result:\n( ( 4 ) ( 5 ) )\n"
                                     "( ( 5 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( sym2, 2UL );
      checkColumns ( sym2, 2UL );
      checkCapacity( sym2, 4UL );
      checkNonZeros( sym2, 4UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );

      if( sym2(0,0) != vec( 1 ) || sym2(0,1) != vec( 2 ) ||
          sym2(1,0) != vec( 2 ) || sym2(1,1) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix swap";

      OST sym1( 2UL );
      sym1(0,0) = vec( 1 );
      sym1(0,1) = vec( 2 );
      sym1(1,1) = vec( 3 );

      OST sym2( 2UL );
      sym2(0,0) = vec( 4 );
      sym2(0,1) = vec( 5 );

      swap( sym1, sym2 );

      checkRows    ( sym1, 2UL );
      checkColumns ( sym1, 2UL );
      checkCapacity( sym1, 4UL );
      checkNonZeros( sym1, 3UL );
      checkNonZeros( sym1, 0UL, 2UL );
      checkNonZeros( sym1, 1UL, 1UL );

      if( sym1(0,0) != vec( 4 ) || sym1(0,1) != vec( 5 ) ||
          sym1(1,0) != vec( 5 ) || !isDefault( sym1(1,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym1 << "\n"
             << "   Expected result:\n( ( 4 ) ( 5 ) )\n"
                                     "( ( 5 ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( sym2, 2UL );
      checkColumns ( sym2, 2UL );
      checkCapacity( sym2, 4UL );
      checkNonZeros( sym2, 4UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );

      if( sym2(0,0) != vec( 1 ) || sym2(0,1) != vec( 2 ) ||
          sym2(1,0) != vec( 2 ) || sym2(1,1) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) )\n"
                                     "( ( 2 ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         ST sym;

         if( isDefault( sym ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         ST sym( 3UL );

         if( isDefault( sym(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sym ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         ST sym( 3UL );
         sym(0,1) = vec( 1 );

         if( isDefault( sym(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sym ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         OST sym;

         if( isDefault( sym ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         OST sym( 3UL );

         if( isDefault( sym(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sym ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         OST sym( 3UL );
         sym(0,1) = vec( 1 );

         if( isDefault( sym(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sym ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
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
void DenseNonNumericTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function (non-overlapping submatrix)";

      typedef blaze::Submatrix<ST>  SMT;

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(2,2) = vec(  3 );

      SMT sm = submatrix( sym, 0UL, 1UL, 2UL, 2UL );

      if( sm(0,1) != vec( 7 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result: " << sm(0,1) << "\n"
             << "   Expected result: ( 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || *it != vec( -4 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: ( -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,1) = vec( -5 );

      if( sm(0,0) != vec( -4 ) || sm(0,1) != vec(  7 ) ||
          sm(1,0) != vec(  2 ) || sm(1,1) != vec( -5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( ( -4 ) (  7 ) )\n"
                                     "( (  2 ) ( -5 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 ) || sym(1,2) != vec( -5 ) ||
          sym(2,0) != vec(  7 ) || sym(2,1) != vec( -5 ) || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) )\n"
                                     "( ( -4 ) (  2 ) ( -5 ) )\n"
                                     "( (  7 ) ( -5 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( !isDefault( sm(0,0) ) || !isDefault( sm(0,1) ) ||
          !isDefault( sm(1,0) ) || !isDefault( sm(1,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( ( ) ( ) )\n( ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || !isDefault( sym(0,2) ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function (non-overlapping submatrix)";

      typedef blaze::Submatrix<OST>  SMT;

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(2,2) = vec(  3 );

      SMT sm = submatrix( sym, 0UL, 1UL, 2UL, 2UL );

      if( sm(0,1) != vec( 7 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result: " << sm(0,1) << "\n"
             << "   Expected result: ( 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || *it != vec( -4 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: ( -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,1) = vec( -5 );

      if( sm(0,0) != vec( -4 ) || sm(0,1) != vec(  7 ) ||
          sm(1,0) != vec(  2 ) || sm(1,1) != vec( -5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( ( -4 ) (  7 ) )\n"
                                     "( (  2 ) ( -5 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 ) || sym(1,2) != vec( -5 ) ||
          sym(2,0) != vec(  7 ) || sym(2,1) != vec( -5 ) || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) )\n"
                                     "( ( -4 ) (  2 ) ( -5 ) )\n"
                                     "( (  7 ) ( -5 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( !isDefault( sm(0,0) ) || !isDefault( sm(0,1) ) ||
          !isDefault( sm(1,0) ) || !isDefault( sm(1,1) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( ( ) ( ) )\n( ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || !isDefault( sym(0,2) ) ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          !isDefault( sym(2,0) ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) (   ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( (   ) (   ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testRow()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      typedef blaze::DenseRow<ST>  RT;

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(2,2) = vec(  3 );

      RT row1 = row( sym, 1UL );

      if( row1[1] != vec( 2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: ( 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it = row1.begin();

      if( it == row1.end() || *it != vec( -4 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: ( -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      row1[2] = vec( -5 );

      if( row1[0] != vec( -4 ) || row1[1] != vec( 2 ) || row1[2] != vec( -5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( ( -4 ) ( 2 ) ( -5 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 ) || sym(1,2) != vec( -5 ) ||
          sym(2,0) != vec(  7 ) || sym(2,1) != vec( -5 ) || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) )\n"
                                     "( ( -4 ) (  2 ) ( -5 ) )\n"
                                     "( (  7 ) ( -5 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( !isDefault( row1[0] ) || !isDefault( row1[1] ) || !isDefault( row1[2] ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( ( ) ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 7 )   ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 7 )   || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 7 ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( ( 7 ) (   ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      typedef blaze::DenseRow<OST>  RT;

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(2,2) = vec(  3 );

      RT row1 = row( sym, 1UL );

      if( row1[1] != vec( 2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: ( 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it = row1.begin();

      if( it == row1.end() || *it != vec( -4 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: ( -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      row1[2] = vec( -5 );

      if( row1[0] != vec( -4 ) || row1[1] != vec( 2 ) || row1[2] != vec( -5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( ( -4 ) ( 2 ) ( -5 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 ) || sym(1,2) != vec( -5 ) ||
          sym(2,0) != vec(  7 ) || sym(2,1) != vec( -5 ) || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) )\n"
                                     "( ( -4 ) (  2 ) ( -5 ) )\n"
                                     "( (  7 ) ( -5 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( !isDefault( row1[0] ) || !isDefault( row1[1] ) || !isDefault( row1[2] ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( ( ) ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 7 )   ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 7 )   || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 7 ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( ( 7 ) (   ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
void DenseNonNumericTest::testColumn()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      typedef blaze::DenseColumn<ST>  CT;

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(2,2) = vec(  3 );

      CT col1 = column( sym, 1UL );

      if( col1[1] != vec( 2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: ( 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it = col1.begin();

      if( it == col1.end() || *it != vec( -4 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: ( -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      col1[2] = vec( -5 );

      if( col1[0] != vec( -4 ) || col1[1] != vec( 2 ) || col1[2] != vec( -5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( ( -4 ) ( 2 ) ( -5 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 ) || sym(1,2) != vec( -5 ) ||
          sym(2,0) != vec(  7 ) || sym(2,1) != vec( -5 ) || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) )\n"
                                     "( ( -4 ) (  2 ) ( -5 ) )\n"
                                     "( (  7 ) ( -5 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col1 );

      if( !isDefault( col1[0] ) || !isDefault( col1[1] ) || !isDefault( col1[2] ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( ( ) ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 7 )   ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 7 )   || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 7 ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( ( 7 ) (   ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      typedef blaze::DenseColumn<OST>  CT;

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(2,2) = vec(  3 );

      CT col1 = column( sym, 1UL );

      if( col1[1] != vec( 2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: ( 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it = col1.begin();

      if( it == col1.end() || *it != vec( -4 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: ( -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      col1[2] = vec( -5 );

      if( col1[0] != vec( -4 ) || col1[1] != vec( 2 ) || col1[2] != vec( -5 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( ( -4 ) ( 2 ) ( -5 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 ) || sym(1,2) != vec( -5 ) ||
          sym(2,0) != vec(  7 ) || sym(2,1) != vec( -5 ) || sym(2,2) != vec(  3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) )\n"
                                     "( ( -4 ) (  2 ) ( -5 ) )\n"
                                     "( (  7 ) ( -5 ) (  3 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col1 );

      if( !isDefault( col1[0] ) || !isDefault( col1[1] ) || !isDefault( col1[2] ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( ( ) ( ) ( ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec( 1 )   || !isDefault( sym(0,1) ) || sym(0,2) != vec( 7 )   ||
          !isDefault( sym(1,0) ) || !isDefault( sym(1,1) ) || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 7 )   || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) (   ) ( 7 ) )\n"
                                     "( (   ) (   ) (   ) )\n"
                                     "( ( 7 ) (   ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
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
   std::cout << "   Running SymmetricMatrix dense non-numeric test..." << std::endl;

   try
   {
      RUN_SYMMETRICMATRIX_DENSENONNUMERIC_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SymmetricMatrix dense non-numeric test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
