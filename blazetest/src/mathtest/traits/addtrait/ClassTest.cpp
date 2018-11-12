//=================================================================================================
/*!
//  \file src/mathtest/traits/addtrait/ClassTest.cpp
//  \brief Source file for the AddTrait class test
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
#include <utility>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/CustomVector.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/HybridVector.h>
#include <blaze/math/IdentityMatrix.h>
#include <blaze/math/InitializerMatrix.h>
#include <blaze/math/InitializerVector.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/StrictlyLowerMatrix.h>
#include <blaze/math/StrictlyUpperMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/UniformMatrix.h>
#include <blaze/math/UniformVector.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/typetraits/Decay.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blazetest/mathtest/traits/addtrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace addtrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the AddTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testScalarAddition();
   testVectorAddition();
   testMatrixAddition();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'AddTrait' class template for scalar/scalar addition operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'AddTrait' class template for scalar/scalar
// addition operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testScalarAddition()
{
   using namespace blaze;


   // int/...
   {
      // .../int
      {
         using T1 = int;
         using T2 = int;
         using RT = int;
         static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );
      }

      // .../double
      {
         using T1 = int;
         using T2 = double;
         using RT = double;
         static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );
      }
   }

   // double/...
   {
      // .../int
      {
         using T1 = double;
         using T2 = int;
         using RT = double;
         static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );
      }

      // .../double
      {
         using T1 = double;
         using T2 = double;
         using RT = double;
         static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );
      }

      // .../complex<double>
      {
         using T1 = double;
         using T2 = complex<double>;
         using RT = complex<double>;
         static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );
      }
   }

   // complex<double>/...
   {
      // .../double
      {
         using T1 = complex<double>;
         using T2 = double;
         using RT = complex<double>;
         static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );
      }

      // .../complex<double>
      {
         using T1 = complex<double>;
         using T2 = complex<double>;
         using RT = complex<double>;
         static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'AddTrait' class template for vector/vector addition operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'AddTrait' class template for vector/vector
// addition operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testVectorAddition()
{
   using namespace blaze;


   // StaticVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = HybridVector<double,5UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = HybridVector<double,5UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CompressedVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = CompressedVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = CompressedVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // HybridVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = HybridVector<double,7UL,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = HybridVector<double,7UL,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CompressedVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = CompressedVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = CompressedVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // DynamicVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = HybridVector<double,7UL,columnVector>;
            using RT = HybridVector<double,7UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = HybridVector<double,7UL,rowVector>;
            using RT = HybridVector<double,7UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CompressedVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = CompressedVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = CompressedVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // CustomVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = HybridVector<double,7UL,columnVector>;
            using RT = HybridVector<double,7UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = HybridVector<double,7UL,rowVector>;
            using RT = HybridVector<double,7UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CompressedVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = CompressedVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = CompressedVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // UniformVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = HybridVector<double,7UL,columnVector>;
            using RT = HybridVector<double,7UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = HybridVector<double,7UL,rowVector>;
            using RT = HybridVector<double,7UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = UniformVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = UniformVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CompressedVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = CompressedVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = CompressedVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // InitializerVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = HybridVector<double,7UL,columnVector>;
            using RT = HybridVector<double,7UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = HybridVector<double,7UL,rowVector>;
            using RT = HybridVector<double,7UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CompressedVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = CompressedVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = CompressedVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // CompressedVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = CompressedVector<int,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CompressedVector<int,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = CompressedVector<int,columnVector>;
            using T2 = HybridVector<double,7UL,columnVector>;
            using RT = HybridVector<double,7UL,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CompressedVector<int,rowVector>;
            using T2 = HybridVector<double,7UL,rowVector>;
            using RT = HybridVector<double,7UL,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = CompressedVector<int,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CompressedVector<int,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = CompressedVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CompressedVector<int,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = CompressedVector<int,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CompressedVector<int,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = CompressedVector<int,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CompressedVector<int,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CompressedVector
      {
         {
            using T1 = CompressedVector<int,columnVector>;
            using T2 = CompressedVector<double,columnVector>;
            using RT = CompressedVector<double,columnVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CompressedVector<int,rowVector>;
            using T2 = CompressedVector<double,rowVector>;
            using RT = CompressedVector<double,rowVector>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'AddTrait' class template for matrix/matrix addition operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'AddTrait' class template for matrix/matrix
// addition operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testMatrixAddition()
{
   using namespace blaze;


   // StaticMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // HybridMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,5UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,5UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,5UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,5UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // DynamicMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // CustomMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UniformMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = UniformMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = UniformMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = UniformMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = UniformMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // InitializerMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // CompressedMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = CompressedMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = CompressedMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CompressedMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // IdentityMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = CompressedMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = CompressedMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DiagonalMatrix< CompressedMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DiagonalMatrix< CompressedMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DiagonalMatrix< CompressedMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DiagonalMatrix< CompressedMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = IdentityMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // SymmetricMatrix/... (real)
   {
      // .../StaticMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // SymmetricMatrix/... (complex)
   {
      // .../StaticMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HybridMatrix<int,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<complex<int>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HybridMatrix<int,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<complex<int>,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HybridMatrix<int,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<complex<int>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HybridMatrix<int,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<complex<int>,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DynamicMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DynamicMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DynamicMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DynamicMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniformMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniformMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniformMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniformMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = InitializerMatrix<int>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = InitializerMatrix<int>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CompressedMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CompressedMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CompressedMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CompressedMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = IdentityMatrix<int,rowMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = IdentityMatrix<int,columnMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = IdentityMatrix<int,rowMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = IdentityMatrix<int,columnMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // HermitianMatrix/... (symmetric)
   {
      // .../StaticMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // HermitianMatrix/... (Hermitian)
   {
      // .../StaticMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HybridMatrix<int,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<complex<int>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HybridMatrix<int,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<complex<int>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HybridMatrix<int,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<complex<int>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HybridMatrix<int,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<complex<int>,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DynamicMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DynamicMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DynamicMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DynamicMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniformMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniformMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniformMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniformMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = InitializerMatrix<int>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = InitializerMatrix<int>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CompressedMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CompressedMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CompressedMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CompressedMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = IdentityMatrix<int,rowMajor>;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = IdentityMatrix<int,columnMajor>;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = IdentityMatrix<int,rowMajor>;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = IdentityMatrix<int,columnMajor>;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // LowerMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UniLowerMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // StrictlyLowerMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UpperMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UniUpperMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // StrictlyUpperMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // DiagonalMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CompressedMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CompressedMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../IdentityMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,rowMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = IdentityMatrix<double,columnMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (real)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix (complex)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (symmetric)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix (Hermitian)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // ... UniLowerMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< AddTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = Decay_t< decltype( std::declval<T1>() + std::declval<T2>() ) >;
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }
}
//*************************************************************************************************

} // namespace addtrait

} // namespace traits

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
   std::cout << "   Running AddTrait class test..." << std::endl;

   try
   {
      RUN_ADDTRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during AddTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
