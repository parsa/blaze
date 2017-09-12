//=================================================================================================
/*!
//  \file blaze/math/views/Forward.h
//  \brief Header file for all forward declarations for views
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_VIEWS_FORWARD_H_
#define _BLAZE_MATH_VIEWS_FORWARD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Forward.h>
#include <blaze/math/views/band/BaseTemplate.h>
#include <blaze/math/views/column/BaseTemplate.h>
#include <blaze/math/views/row/BaseTemplate.h>
#include <blaze/math/views/submatrix/BaseTemplate.h>
#include <blaze/math/views/subvector/BaseTemplate.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< bool AF, size_t I, size_t N, typename VT, bool TF >
decltype(auto) subvector( Vector<VT,TF>& );

template< bool AF, size_t I, size_t N, typename VT, bool TF >
decltype(auto) subvector( const Vector<VT,TF>& );

template< bool AF, size_t I, size_t N, typename VT, bool TF >
decltype(auto) subvector( Vector<VT,TF>&& );

template< bool AF, typename VT, bool TF >
decltype(auto) subvector( Vector<VT,TF>&, size_t, size_t );

template< bool AF, typename VT, bool TF >
decltype(auto) subvector( const Vector<VT,TF>&, size_t, size_t );

template< bool AF, typename VT, bool TF >
decltype(auto) subvector( Vector<VT,TF>&&, size_t, size_t );

template< bool AF, size_t I, size_t J, size_t M, size_t N, typename MT, bool SO >
decltype(auto) submatrix( Matrix<MT,SO>& );

template< bool AF, size_t I, size_t J, size_t M, size_t N, typename MT, bool SO >
decltype(auto) submatrix( const Matrix<MT,SO>& );

template< bool AF, size_t I, size_t J, size_t M, size_t N, typename MT, bool SO >
decltype(auto) submatrix( Matrix<MT,SO>&& );

template< bool AF, typename MT, bool SO >
decltype(auto) submatrix( Matrix<MT,SO>&, size_t, size_t, size_t, size_t );

template< bool AF, typename MT, bool SO >
decltype(auto) submatrix( const Matrix<MT,SO>&, size_t, size_t, size_t, size_t );

template< bool AF, typename MT, bool SO >
decltype(auto) submatrix( Matrix<MT,SO>&&, size_t, size_t, size_t, size_t );

template< size_t I, typename MT, bool SO >
decltype(auto) row( Matrix<MT,SO>& );

template< size_t I, typename MT, bool SO >
decltype(auto) row( const Matrix<MT,SO>& );

template< size_t I, typename MT, bool SO >
decltype(auto) row( Matrix<MT,SO>&& );

template< typename MT, bool SO >
decltype(auto) row( Matrix<MT,SO>&, size_t );

template< typename MT, bool SO >
decltype(auto) row( const Matrix<MT,SO>&, size_t );

template< typename MT, bool SO >
decltype(auto) row( Matrix<MT,SO>&&, size_t );

template< size_t I, typename MT, bool SO >
Column<MT,I> column( Matrix<MT,SO>& );

template< size_t I, typename MT, bool SO >
const Column<const MT,I> column( const Matrix<MT,SO>& );

template< size_t I, typename MT, bool SO >
Column<MT,I> column( Matrix<MT,SO>&& );

template< typename MT, bool SO >
Column<MT> column( Matrix<MT,SO>&, size_t );

template< typename MT, bool SO >
const Column<const MT> column( const Matrix<MT,SO>&, size_t );

template< typename MT, bool SO >
Column<MT> column( Matrix<MT,SO>&&, size_t );

template< ptrdiff_t I, typename MT, bool SO >
Band<MT,I> band( Matrix<MT,SO>& );

template< ptrdiff_t I, typename MT, bool SO >
const Band<const MT,I> band( const Matrix<MT,SO>& );

template< ptrdiff_t I, typename MT, bool SO >
Band<MT,I> band( Matrix<MT,SO>&& );

template< typename MT, bool SO >
Band<MT> band( Matrix<MT,SO>&, ptrdiff_t );

template< typename MT, bool SO >
const Band<const MT> band( const Matrix<MT,SO>&, ptrdiff_t );

template< typename MT, bool SO >
Band<MT> band( Matrix<MT,SO>&&, ptrdiff_t );

template< typename MT, bool SO >
decltype(auto) diagonal( Matrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) diagonal( const Matrix<MT,SO>& );

template< typename MT, bool SO >
decltype(auto) diagonal( Matrix<MT,SO>&& );

} // namespace blaze

#endif
