//=================================================================================================
/*!
//  \file blazetest/mathtest/RemoveIdentity.h
//  \brief Header file for the RemoveIdentity type trait
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

#ifndef _BLAZETEST_MATHTEST_REMOVEIDENTITY_H_
#define _BLAZETEST_MATHTEST_REMOVEIDENTITY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/adaptors/Forward.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/util/mpl/If.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Removes the identity property of the given matrix type.
//
// The RemoveIdentity type trait removes the identity property of the given matrix type \a T. In
// case \a T is an identity matrix, the resulting type \a Type will be a compressed DiagonalMatrix
// with the according element type. Else \a Type is set to the given type \a T. Note that this type
// trait only works for matrix types. The attempt to instantiate it with non-matrix types results
// in a compile time error.
*/
template< typename T >  // The given matrix type
struct RemoveIdentity
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using ET = blaze::ElementType_t<T>;
   static constexpr bool SO = blaze::IsColumnMajorMatrix_v<T>;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = blaze::If_t< blaze::IsIdentity_v<T>
                           , blaze::DiagonalMatrix< blaze::CompressedMatrix<ET,SO> >
                           , T >;
   /*! \endcond */
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( T );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the RemoveIdentity class template.
//
// The RemoveIdentity_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the RemoveIdentity class template. For instance, given the types \a T1 and \a T2
// the following two type definitions are identical:

   \code
   using Type1 = typename RemoveIdentity<T>::Type;
   using Type2 = RemoveIdentity_t<T>;
   \endcode
*/
template< typename T >  // The given matrix type
using RemoveIdentity_t = typename RemoveIdentity<T>::Type;
//*************************************************************************************************

} // namespace blazetest

#endif
