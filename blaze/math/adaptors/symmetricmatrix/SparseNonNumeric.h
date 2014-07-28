//=================================================================================================
/*!
//  \file blaze/math/adaptors/symmetricmatrix/SparseNonNumeric.h
//  \brief SymmetricMatrix specialization for sparse matrices with non-numeric element type
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

#ifndef _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_SPARSENONNUMERIC_H_
#define _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_SPARSENONNUMERIC_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/adaptors/symmetricmatrix/BaseTemplate.h>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/Resizable.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SPARSE MATRICES WITH NON-NUMERIC ELEMENT TYPE
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of SymmetricMatrix for sparse matrices with non-numeric element type.
// \ingroup symmetric_matrix
//
// This specialization of SymmetricMatrix adapts the class template to the requirements of sparse
// matrices with non-numeric element type.
*/
template< typename MT >  // Type of the adapted sparse matrix
class SymmetricMatrix<MT,false,false>
   : public SparseMatrix< SymmetricMatrix<MT,false,false>, IsColumnMajorMatrix<MT>::value >
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT::OppositeType   OT;  //!< Opposite type of the sparse matrix.
   typedef typename MT::TransposeType  TT;  //!< Transpose type of the sparse matrix.
   typedef typename MT::ElementType    ET;  //!< Element type of the sparse matrix.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SymmetricMatrix<MT,false,false>  This;            //!< Type of this SymmetricMatrix instance.
   typedef This                             ResultType;      //!< Result type for expression template evaluations.
   typedef SymmetricMatrix<OT,false,false>  OppositeType;    //!< Result type with opposite storage order for expression template evaluations.
   typedef SymmetricMatrix<TT,false,false>  TransposeType;   //!< Transpose type for expression template evaluations.
   typedef ET                               ElementType;     //!< Type of the matrix elements.
   typedef typename MT::ReturnType          ReturnType;      //!< Return type for expression template evaluations.
   typedef const This&                      CompositeType;   //!< Data type for composite expression templates.
   typedef SymmetricProxy<MT>               Reference;       //!< Reference to a non-constant matrix value.
   typedef typename MT::ConstReference      ConstReference;  //!< Reference to a constant matrix value.
   typedef typename MT::ConstIterator       ConstIterator;   //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST          ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_EXPRESSION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_NUMERIC_TYPE   ( ElementType );
   BLAZE_STATIC_ASSERT( IsResizable<MT>::value || IsSquare<MT>::value );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
