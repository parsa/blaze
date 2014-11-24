//=================================================================================================
/*!
//  \file blaze/math/proxy/Proxy.h
//  \brief Header file for the Proxy class
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

#ifndef _BLAZE_MATH_PROXY_PROXY_H_
#define _BLAZE_MATH_PROXY_PROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/proxy/ComplexProxy.h>
#include <blaze/math/proxy/DefaultProxy.h>
#include <blaze/math/proxy/DenseMatrixProxy.h>
#include <blaze/math/proxy/DenseVectorProxy.h>
#include <blaze/math/proxy/SparseMatrixProxy.h>
#include <blaze/math/proxy/SparseVectorProxy.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Proxy base class.
// \ingroup math
//
// The Proxy class is a base class for all proxy classes within the \b Blaze library that may
// represent non-numeric data types (vectors, matrices, ...). It augments the interface of the
// deriving proxy class depending on the data type represented by the proxy. In addition, it
// provides an abstraction from the actual type of the proxy, but enables a type-safe conversion
// back to this type via the 'Curiously Recurring Template Pattern' (CRTP).
//
// In order to use the Proxy class it is necessary to publicly derive from it and to provide
// an accessible member function called \a get(), which grants access to the represented element
// via non-const reference. The following example demonstrates these requirements by means of
// the VectorAccessProxy class:

   \code
   template< typename VT >
   class VectorAccessProxy : public Proxy< VectorAccessProxy<VT>, typename VT::ElementType >
   {
      // ...
      typedef typename VT::ElementType  RepresentedType;
      inline RepresentedType& get() const;
      // ...
   };
   \endcode

// The first template parameter specifies the type of the deriving proxy class (CRTP), the second
// template parameter specifies the type of the element represented by the proxy. Within the
// context of the VectorAccessProxy this is the type of the elements of the vector to be accessed.
// Depending on this type the proxy selects the additional interface to provide to the deriving
// class.
*/
template< typename PT    // Type of the proxy
        , typename RT >  // Type of the represented element
class Proxy : public If< IsVector<RT>
                       , typename If< IsDenseVector<RT>
                                    , DenseVectorProxy<PT,RT>
                                    , SparseVectorProxy<PT,RT>
                                    >::Type
                       , typename If< IsMatrix<RT>
                                    , typename If< IsDenseMatrix<RT>
                                                 , DenseMatrixProxy<PT,RT>
                                                 , SparseMatrixProxy<PT,RT>
                                                 >::Type
                                    , typename If< IsComplex<RT>
                                                 , ComplexProxy<PT,RT>
                                                 , DefaultProxy<PT,RT>
                                                 >::Type
                                    >::Type
                       >::Type
{};
//*************************************************************************************************

} // namespace blaze

#endif
