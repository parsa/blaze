//=================================================================================================
/*!
//  \file blaze/math/typetraits/UnrestrictedType.h
//  \brief Header file for the UnrestrictedType type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_UNRESTRICTEDTYPE_H_
#define _BLAZE_MATH_TYPETRAITS_UNRESTRICTEDTYPE_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the associated type with unrestricted data access.
// \ingroup math_type_traits
//
// Via this type trait it is possible to determine the associated data type with unrestricted data
// access of the given type \a T. In case the given data type has a restriction on its data access
// the returned type provides the same interface, but is free of any restrictions. In case the
// given data type does not have any restrictions on its data access, the returned type is \a T.
// Note that cv-qualifiers are preserved:

   \code
   typedef blaze::LowerMatrix< blaze::DynamicMatrix<double> >  Lower;
   typedef blaze::LowerMatrix< blaze::CompressedMatrix<int> >  Upper;

   blaze::UnrestrictedType< Lower >::Type           // Results in 'blaze::DynamicMatrix<double>'
   blaze::UnrestrictedType< const Upper >::Type     // Results in 'const blaze::CompressedMatrix<int>'
   blaze::UnrestrictedType< volatile Lower >::Type  // Results in 'volatile blaze::DynamicMatrix<double>'
   blaze::UnrestrictedType< int >::Type             // Results in 'int'
   \endcode
*/
template< typename T >
struct UnrestrictedType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef T  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UnrestrictedType type trait for const types.
// \ingroup math_type_traits
*/
template< typename T >
struct UnrestrictedType< const T >
{
 public:
   //**********************************************************************************************
   typedef const typename UnrestrictedType<T>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UnrestrictedType type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename T >
struct UnrestrictedType< volatile T >
{
 public:
   //**********************************************************************************************
   typedef volatile typename UnrestrictedType<T>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UnrestrictedType type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename T >
struct UnrestrictedType< const volatile T >
{
 public:
   //**********************************************************************************************
   typedef const volatile typename UnrestrictedType<T>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
