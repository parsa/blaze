//=================================================================================================
/*!
//  \file blaze/math/functors/Bind1st.h
//  \brief Header file for the Bind1st functor
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_FUNCTORS_BIND1ST_H_
#define _BLAZE_MATH_FUNCTORS_BIND1ST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/SIMDPack.h>
#include <blaze/math/simd/Set.h>
#include <blaze/math/typetraits/IsSIMDEnabled.h>
#include <blaze/math/typetraits/YieldsSymmetric.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/system/Inline.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for a binary operation with fixed 1st argument.
// \ingroup functors
*/
template< typename OP    // Type of the binary operation
        , typename A1 >  // Type of the bound argument
struct Bind1st
{
 public:
   //**********************************************************************************************
   /*!\brief Constructor of the Bind1st functor.
   //
   // \param op The binary operation.
   // \param a1 The 1st argument.
   */
   explicit inline constexpr Bind1st( const OP& op, const A1& a1 )
      : op_( op )  // The wrapped operation.
      , a1_( a1 )  // The 1st argument
   {}
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given object/value.
   //
   // \param a The given object/value.
   // \return The result of the wrapped operation for the given object/value.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE constexpr decltype(auto) operator()( const T& a ) const
   {
      return op_( a1_, a );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data type \a T.
   //
   // \return \a true in case SIMD is enabled for the data type \a T, \a false if not.
   */
   template< typename T >
   static constexpr bool simdEnabled() { return IsSIMDEnabled_v<OP,A1,T>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operation supports padding, i.e. whether it can deal with zeros.
   //
   // \return \a true in case padding is supported, \a false if not.
   */
   static constexpr bool paddingEnabled() { return false; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given SIMD vector.
   //
   // \param a The given SIMD vector.
   // \return The result of the wrapped operation for the given SIMD vector.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const T& a ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T );
      return op_.load( set( a1_ ), a );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   OP op_;  //!< The binary operation.
   A1 a1_;  //!< The 1st argument.
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Binds the given object/value to the 1st parameter of the given binary operation.
// \ingroup functors
//
// \param op The binary operation to be wrapped.
// \param x The argument to be bound to the second parameter of the binary operation.
// \return The operation with bound 1st argument.
//
// The \a bind1st() function binds the given argument \a x to the 1st parameter of the given
// binary operation \a op.
*/
template< typename OP    // Type of the binary operation
        , typename A1 >  // Type of the bound argument
inline constexpr Bind1st<OP,A1> bind1st( const OP& op, const A1& a1 )
{
   return Bind1st<OP,A1>( op, a1 );
}
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename A1, typename T >
struct YieldsUniform<Bind1st<OP,A1>,T>
   : public YieldsUniform<OP,T>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename A1, typename MT >
struct YieldsSymmetric<Bind1st<OP,A1>,MT>
   : public YieldsSymmetric<OP,MT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
