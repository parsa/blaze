//=================================================================================================
/*!
//  \file blaze/math/traits/DeclLowTrait.h
//  \brief Header file for the decllow trait
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

#ifndef _BLAZE_MATH_TRAITS_DECLLOWTRAIT_H_
#define _BLAZE_MATH_TRAITS_DECLLOWTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/adaptors/lowermatrix/BaseTemplate.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/Decay.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the DeclLowTrait class.
// \ingroup math_traits
//
// \section decllowtrait_general General
//
// The DeclLowTrait class template offers the possibility to select the resulting data type
// of a generic decllow() operation on the given type \a MT. DeclLowTrait defines the nested
// type \a Type, which represents the resulting data type of the decllow() operation. In case
// the given data type is not a dense or sparse matrix type, the resulting data type \a Type is
// set to \a INVALID_TYPE. Note that \a const and \a volatile qualifiers and reference modifiers
// are generally ignored.
//
// Per default, the DeclLowTrait template only supports the following matrix types:
//
// <ul>
//    <li>blaze::StaticMatrix</li>
//    <li>blaze::HybridMatrix</li>
//    <li>blaze::DynamicMatrix</li>
//    <li>blaze::CustomMatrix</li>
//    <li>blaze::CompressedMatrix</li>
//    <li>blaze::IdentityMatrix</li>
//    <li>blaze::SymmetricMatrix</li>
//    <li>blaze::HermitianMatrix</li>
//    <li>blaze::LowerMatrix</li>
//    <li>blaze::UniLowerMatrix</li>
//    <li>blaze::StrictlyLowerMatrix</li>
//    <li>blaze::UpperMatrix</li>
//    <li>blaze::UniUpperMatrix</li>
//    <li>blaze::StrictlyUpperMatrix</li>
//    <li>blaze::DiagonalMatrix</li>
//    <li>blaze::Submatrix</li>
// </ul>
//
//
// \section decllowtrait_specializations Creating custom specializations
//
// It is possible to specialize the DeclLowTrait template for additional user-defined matrix types.
// The following example shows the according specialization for the SymmetricMatrix class template:

   \code
   template< typename MT, bool SO, bool DF, bool NF >
   struct DeclLowTrait< SymmetricMatrix<MT,SO,DF,NF> >
   {
      using Type = DiagonalMatrix<MT>;
   };
   \endcode

// \n \section decllowtrait_examples Examples
//
// The following example demonstrates the use of the DeclLowTrait template, where depending on
// the given matrix type the resulting type is selected:

   \code
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;
   using blaze::DeclLowTrait;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the resulting type of a row-major dynamic matrix
   using MatrixType1  = DynamicMatrix<int,rowMajor>;
   using DeclLowType1 = typename DeclLowTrait<MatrixType1>::Type;

   // Definition of the resulting type of a symmetric column-major static matrix
   using MatrixType2  = SymmetricMatrix< StaticMatrix<int,3UL,3UL,columnMajor> >;
   using DeclLowType2 = typename DeclLowTrait<MatrixType2>::Type;
   \endcode
*/
template< typename MT >  // Type of the matrix
struct DeclLowTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { using Type = INVALID_TYPE; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Result { using Type = LowerMatrix<MT>; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<MT>, IsVolatile<MT>, IsReference<MT> >
                            , DeclLowTrait< Decay_<MT> >
                            , If_< IsMatrix<MT>
                                 , Result
                                 , Failure > >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the DeclLowTrait type trait.
// \ingroup math_traits
//
// The DeclLowTrait_ alias declaration provides a convenient shortcut to access the nested
// \a Type of the DeclLowTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename DeclLowTrait<MT>::Type;
   using Type2 = DeclLowTrait_<MT>;
   \endcode
*/
template< typename MT >  // Type of the matrix
using DeclLowTrait_ = typename DeclLowTrait<MT>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
