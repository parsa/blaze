//=================================================================================================
/*!
//  \file src/mathtest/typetraits/OperationTest.cpp
//  \brief Source file for the mathematical type traits operation test
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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/Diagonal.h>
#include <blaze/math/constraints/Identity.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/UniLower.h>
#include <blaze/math/constraints/UniUpper.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/constraints/Vector.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingNumeric.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/mathtest/typetraits/OperationTest.h>


namespace blazetest {

namespace mathtest {

namespace typetraits {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the OperationTest class test.
//
// \exception std::runtime_error Operation error detected.
*/
OperationTest::OperationTest()
{
   testIsColumnMajorMatrix();
   testIsColumnVector();
   testIsDiagonal();
   testIsIdentity();
   testIsLower();
   testIsMatrix();
   testIsRowVector();
   testIsSymmetric();
   testIsUniLower();
   testIsUniUpper();
   testIsUpper();
   testIsVector();
   testRemoveAdaptor();
   testUnderlyingBuiltin();
   testUnderlyingNumeric();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST TYPE TRAITS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the mathematical 'testIsColumnMajorMatrix' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsColumnMajorMatrix' type
// trait. In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsColumnMajorMatrix()
{
   using blaze::rowMajor;
   using blaze::columnMajor;

   typedef blaze::StaticMatrix<float,3U,3U,columnMajor>       Type1;
   typedef const blaze::DynamicMatrix<double,columnMajor>     Type2;
   typedef volatile blaze::CompressedMatrix<int,columnMajor>  Type3;
   typedef blaze::StaticMatrix<float,3U,3U,rowMajor>          Type4;
   typedef const blaze::DynamicMatrix<double,rowMajor>        Type5;
   typedef volatile blaze::CompressedMatrix<int,rowMajor>     Type6;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type2 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsColumnVector' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsColumnVector' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsColumnVector()
{
   using blaze::columnVector;
   using blaze::rowVector;

   typedef blaze::StaticVector<float,3U,columnVector>          Type1;
   typedef const blaze::DynamicVector<double,columnVector>     Type2;
   typedef volatile blaze::CompressedVector<int,columnVector>  Type3;
   typedef blaze::StaticVector<float,3U,rowVector>             Type4;
   typedef const blaze::DynamicVector<double,rowVector>        Type5;
   typedef volatile blaze::CompressedVector<int,rowVector>     Type6;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type2 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsDiagonal' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsDiagonal' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsDiagonal()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;

   typedef SymmetricMatrix< DynamicMatrix<int> >           Type1;
   typedef const SymmetricMatrix< DynamicMatrix<int> >     Type2;
   typedef volatile SymmetricMatrix< DynamicMatrix<int> >  Type3;
   typedef LowerMatrix< DynamicMatrix<int> >               Type4;
   typedef const LowerMatrix< DynamicMatrix<int> >         Type5;
   typedef volatile LowerMatrix< DynamicMatrix<int> >      Type6;
   typedef UpperMatrix< DynamicMatrix<int> >               Type7;
   typedef const UpperMatrix< DynamicMatrix<int> >         Type8;
   typedef volatile UpperMatrix< DynamicMatrix<int> >      Type9;
   typedef DiagonalMatrix< DynamicMatrix<int> >            Type10;
   typedef const DiagonalMatrix< DynamicMatrix<int> >      Type11;
   typedef volatile DiagonalMatrix< DynamicMatrix<int> >   Type12;

   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type4  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type5  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type6  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type7  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type8  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type9  );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type10 );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type11 );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type12 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsIdentity' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsIdentity' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsIdentity()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;

   typedef SymmetricMatrix< DynamicMatrix<int> >           Type1;
   typedef const SymmetricMatrix< DynamicMatrix<int> >     Type2;
   typedef volatile SymmetricMatrix< DynamicMatrix<int> >  Type3;
   typedef LowerMatrix< DynamicMatrix<int> >               Type4;
   typedef const LowerMatrix< DynamicMatrix<int> >         Type5;
   typedef volatile LowerMatrix< DynamicMatrix<int> >      Type6;
   typedef UpperMatrix< DynamicMatrix<int> >               Type7;
   typedef const UpperMatrix< DynamicMatrix<int> >         Type8;
   typedef volatile UpperMatrix< DynamicMatrix<int> >      Type9;
   typedef DiagonalMatrix< DynamicMatrix<int> >            Type10;
   typedef const DiagonalMatrix< DynamicMatrix<int> >      Type11;
   typedef volatile DiagonalMatrix< DynamicMatrix<int> >   Type12;

   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type5  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type6  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type7  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type8  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type9  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type10 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type11 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type12 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsLower' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsLower' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsLower()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;

   typedef SymmetricMatrix< DynamicMatrix<int> >           Type1;
   typedef const SymmetricMatrix< DynamicMatrix<int> >     Type2;
   typedef volatile SymmetricMatrix< DynamicMatrix<int> >  Type3;
   typedef LowerMatrix< DynamicMatrix<int> >               Type4;
   typedef const LowerMatrix< DynamicMatrix<int> >         Type5;
   typedef volatile LowerMatrix< DynamicMatrix<int> >      Type6;
   typedef UpperMatrix< DynamicMatrix<int> >               Type7;
   typedef const UpperMatrix< DynamicMatrix<int> >         Type8;
   typedef volatile UpperMatrix< DynamicMatrix<int> >      Type9;
   typedef DiagonalMatrix< DynamicMatrix<int> >            Type10;
   typedef const DiagonalMatrix< DynamicMatrix<int> >      Type11;
   typedef volatile DiagonalMatrix< DynamicMatrix<int> >   Type12;

   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type3  );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type4  );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type5  );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type6  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type7  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type8  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type9  );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type10 );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type11 );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type12 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsMatrix' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsMatrix' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsMatrix()
{
   typedef blaze::StaticMatrix<float,3U,3U,false>      Type1;
   typedef const blaze::DynamicMatrix<double,true>     Type2;
   typedef volatile blaze::CompressedMatrix<int,true>  Type3;
   typedef blaze::StaticVector<float,3U,false>         Type4;
   typedef const blaze::DynamicVector<double,true>     Type5;
   typedef volatile blaze::CompressedVector<int,true>  Type6;

   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type1 );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type2 );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsRowMajorMatrix' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsRowMajorMatrix' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsRowMajorMatrix()
{
   using blaze::rowMajor;
   using blaze::columnMajor;

   typedef blaze::StaticMatrix<float,3U,3U,rowMajor>          Type1;
   typedef const blaze::DynamicMatrix<double,rowMajor>        Type2;
   typedef volatile blaze::CompressedMatrix<int,rowMajor>     Type3;
   typedef blaze::StaticMatrix<float,3U,3U,columnMajor>       Type4;
   typedef const blaze::DynamicMatrix<double,columnMajor>     Type5;
   typedef volatile blaze::CompressedMatrix<int,columnMajor>  Type6;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsRowVector' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsRowVector' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsRowVector()
{
   using blaze::columnVector;
   using blaze::rowVector;

   typedef blaze::StaticVector<float,3U,rowVector>             Type1;
   typedef const blaze::DynamicVector<double,rowVector>        Type2;
   typedef volatile blaze::CompressedVector<int,rowVector>     Type3;
   typedef blaze::StaticVector<float,3U,columnVector>          Type4;
   typedef const blaze::DynamicVector<double,columnVector>     Type5;
   typedef volatile blaze::CompressedVector<int,columnVector>  Type6;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsSymmetric' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsSymmetric' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsSymmetric()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;

   typedef SymmetricMatrix< DynamicMatrix<int> >           Type1;
   typedef const SymmetricMatrix< DynamicMatrix<int> >     Type2;
   typedef volatile SymmetricMatrix< DynamicMatrix<int> >  Type3;
   typedef LowerMatrix< DynamicMatrix<int> >               Type4;
   typedef const LowerMatrix< DynamicMatrix<int> >         Type5;
   typedef volatile LowerMatrix< DynamicMatrix<int> >      Type6;
   typedef UpperMatrix< DynamicMatrix<int> >               Type7;
   typedef const UpperMatrix< DynamicMatrix<int> >         Type8;
   typedef volatile UpperMatrix< DynamicMatrix<int> >      Type9;
   typedef DiagonalMatrix< DynamicMatrix<int> >            Type10;
   typedef const DiagonalMatrix< DynamicMatrix<int> >      Type11;
   typedef volatile DiagonalMatrix< DynamicMatrix<int> >   Type12;

   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type1  );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type2  );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type3  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type5  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type6  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type7  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type8  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type9  );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type10 );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type11 );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type12 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsUniLower' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsUniLower' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsUniLower()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;

   typedef SymmetricMatrix< DynamicMatrix<int> >           Type1;
   typedef const SymmetricMatrix< DynamicMatrix<int> >     Type2;
   typedef volatile SymmetricMatrix< DynamicMatrix<int> >  Type3;
   typedef LowerMatrix< DynamicMatrix<int> >               Type4;
   typedef const LowerMatrix< DynamicMatrix<int> >         Type5;
   typedef volatile LowerMatrix< DynamicMatrix<int> >      Type6;
   typedef UpperMatrix< DynamicMatrix<int> >               Type7;
   typedef const UpperMatrix< DynamicMatrix<int> >         Type8;
   typedef volatile UpperMatrix< DynamicMatrix<int> >      Type9;
   typedef DiagonalMatrix< DynamicMatrix<int> >            Type10;
   typedef const DiagonalMatrix< DynamicMatrix<int> >      Type11;
   typedef volatile DiagonalMatrix< DynamicMatrix<int> >   Type12;

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type3  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type6  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type7  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type8  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type9  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type10 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type11 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type12 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsUniUpper' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsUniUpper' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsUniUpper()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;

   typedef SymmetricMatrix< DynamicMatrix<int> >           Type1;
   typedef const SymmetricMatrix< DynamicMatrix<int> >     Type2;
   typedef volatile SymmetricMatrix< DynamicMatrix<int> >  Type3;
   typedef LowerMatrix< DynamicMatrix<int> >               Type4;
   typedef const LowerMatrix< DynamicMatrix<int> >         Type5;
   typedef volatile LowerMatrix< DynamicMatrix<int> >      Type6;
   typedef UpperMatrix< DynamicMatrix<int> >               Type7;
   typedef const UpperMatrix< DynamicMatrix<int> >         Type8;
   typedef volatile UpperMatrix< DynamicMatrix<int> >      Type9;
   typedef DiagonalMatrix< DynamicMatrix<int> >            Type10;
   typedef const DiagonalMatrix< DynamicMatrix<int> >      Type11;
   typedef volatile DiagonalMatrix< DynamicMatrix<int> >   Type12;

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type4  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type6  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type7  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type8  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type9  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type10 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type11 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type12 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsUpper' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsUpper' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsUpper()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;

   typedef SymmetricMatrix< DynamicMatrix<int> >           Type1;
   typedef const SymmetricMatrix< DynamicMatrix<int> >     Type2;
   typedef volatile SymmetricMatrix< DynamicMatrix<int> >  Type3;
   typedef LowerMatrix< DynamicMatrix<int> >               Type4;
   typedef const LowerMatrix< DynamicMatrix<int> >         Type5;
   typedef volatile LowerMatrix< DynamicMatrix<int> >      Type6;
   typedef UpperMatrix< DynamicMatrix<int> >               Type7;
   typedef const UpperMatrix< DynamicMatrix<int> >         Type8;
   typedef volatile UpperMatrix< DynamicMatrix<int> >      Type9;
   typedef DiagonalMatrix< DynamicMatrix<int> >            Type10;
   typedef const DiagonalMatrix< DynamicMatrix<int> >      Type11;
   typedef volatile DiagonalMatrix< DynamicMatrix<int> >   Type12;

   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type4  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type5  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type6  );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type7  );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type8  );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type9  );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type10 );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type11 );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type12 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsVector' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsVector' type trait. In case
// an error is detected, a compilation error is created.
*/
void OperationTest::testIsVector()
{
   typedef blaze::StaticVector<float,3U,false>         Type1;
   typedef const blaze::DynamicVector<double,true>     Type2;
   typedef volatile blaze::CompressedVector<int,true>  Type3;
   typedef blaze::StaticMatrix<double,3U,3U,false>     Type4;
   typedef const blaze::DynamicMatrix<double,true>     Type5;
   typedef volatile blaze::CompressedMatrix<int,true>  Type6;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type1 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type2 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'RemoveAdaptor' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'RemoveAdaptor' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testRemoveAdaptor()
{
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;
   using blaze::SymmetricMatrix;
   using blaze::RemoveAdaptor;

   typedef SymmetricMatrix< DynamicMatrix<int> >              Source1;
   typedef const SymmetricMatrix< CompressedMatrix<float> >   Source2;
   typedef volatile SymmetricMatrix< DynamicMatrix<double> >  Source3;
   typedef int                                                Source4;
   typedef const DynamicVector<int>                           Source5;
   typedef volatile DynamicMatrix<int>                        Source6;

   typedef DynamicMatrix<int>              Result1;
   typedef const CompressedMatrix<float>   Result2;
   typedef volatile DynamicMatrix<double>  Result3;
   typedef int                             Result4;
   typedef const DynamicVector<int>        Result5;
   typedef volatile DynamicMatrix<int>     Result6;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source1>::Type, Result1 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source2>::Type, Result2 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source3>::Type, Result3 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source4>::Type, Result4 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source5>::Type, Result5 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source6>::Type, Result6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'UnderlyingBuiltin' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'UnderlyingBuiltin' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testUnderlyingBuiltin()
{
   using blaze::complex;
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::UnderlyingBuiltin;

   typedef double                                    Type1;  // Built-in data type
   typedef complex<float>                            Type2;  // Complex data type
   typedef StaticVector<int,3UL>                     Type3;  // Vector with built-in element type
   typedef CompressedVector< DynamicVector<float> >  Type4;  // Vector with vector element type

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1>::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2>::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4>::Type, float );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'UnderlyingNumeric' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'UnderlyingNumeric' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testUnderlyingNumeric()
{
   using blaze::complex;
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::UnderlyingNumeric;

   typedef double                                    Type1;  // Built-in data type
   typedef complex<float>                            Type2;  // Complex data type
   typedef StaticVector<int,3UL>                     Type3;  // Vector with built-in element type
   typedef CompressedVector< DynamicVector<float> >  Type4;  // Vector with vector element type

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1>::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2>::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4>::Type, float );
}
//*************************************************************************************************

} // namespace typetraits

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
   std::cout << "   Running mathematical type traits operation test..." << std::endl;

   try
   {
      RUN_TYPETRAITS_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during mathematical type traits operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
