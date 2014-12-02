//=================================================================================================
/*!
//  \file blaze/util/Null.h
//  \brief Header file for a safe C++ NULL pointer implementation
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

#ifndef _BLAZE_UTIL_NULL_H_
#define _BLAZE_UTIL_NULL_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Null.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Safe C++ NULL pointer implementation.
// \ingroup util
//
// This implementation offers a remedy for the use of the NULL pointer in C++. For this, the
// NULL macro is replaced by an instance of the Null class, which can only be assigned and
// compared with pointers and pointers-to-member. Therefore the use of NULL regains the type
// safety it lost in C++ due to the strict C++ type system.\n
// The NULL pointer is used exactly as before:

   \code
   int* pi = NULL;
   if( pi == NULL ) {...}
   \endcode
*/
class Null
{
 public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   inline Null();
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   template< typename T >
   inline operator T*() const;

   template< typename T, typename C >
   inline operator T C::*() const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename T >
   inline bool equal( const T* rhs ) const;

   template< typename T, typename C >
   inline bool equal( const T C::* rhs ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   Null( const Null& n );             //!< Copy constructor (private & undefined)
   Null& operator=( const Null& n );  //!< Copy assignment operator (private & undefined)
   void* operator&() const;           //!< Address operator (private & undefined)
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor of the Null class.
*/
inline Null::Null()
{}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion operator to a pointer.
//
// This conversion operator offers a type safe conversion of zero to a pointer of any kind.
*/
template< typename T >
inline Null::operator T*() const
{
   return 0;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to a pointer-to-member.
//
// This conversion operator offers the type safe conversion of zero to a pointer-to-member of
// any kind.
*/
template< typename T, typename C >
inline Null::operator T C::*() const
{
   return 0;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Comparison between Null and a pointer.
//
// The function offers a type safe comparison between zero and an arbitrary pointer.
*/
template< typename T >
inline bool Null::equal( const T* rhs ) const
{
   return rhs == 0;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Comparison between Null and a pointer-to-member.
//
// The function offers a type safe comparison between zero and an arbitrary pointer-to-member.
*/
template< typename T, typename C >
inline bool Null::equal( const T C::* rhs ) const
{
   return rhs == 0;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Null operators */
//@{
template< typename T >
inline bool operator==( const Null& lhs, const T& rhs );

template< typename T >
inline bool operator==( const T& lhs, const Null& rhs );

template< typename T >
inline bool operator!=( const Null& lhs, const T& rhs );

template< typename T >
inline bool operator!=( const T& lhs, const Null& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between Null and a pointer or pointer-to-member.
//
// This operator takes a reference to an object of type T instead of a pointer of a pointer-
// to-member to avoid the ambiguity with the built-in pointer comparison operators. However,
// only pointers and pointers-to-member can be compared to Null.
*/
template< typename T >
inline bool operator==( const Null& lhs, const T& rhs )
{
   return lhs.equal( rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a pointer or pointer-to-member and Null.
//
// This operator takes a reference to an object of type T instead of a pointer of a pointer-
// to-member to avoid the ambiguity with the built-in pointer comparison operators. However,
// only pointers and pointers-to-member can be compared to Null.
*/
template< typename T >
inline bool operator==( const T& lhs, const Null& rhs )
{
   return rhs.equal( lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between Null and a pointer or pointer-to-member.
//
// This operator takes a reference to an object of type T instead of a pointer of a pointer-
// to-member to avoid the ambiguity with the built-in pointer comparison operators. However,
// only pointers and pointers-to-member can be compared to Null.
*/
template< typename T >
inline bool operator!=( const Null& lhs, const T& rhs )
{
   return !lhs.equal( rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a pointer or pointer-to-member and Null.
//
// This operator takes a reference to an object of type T instead of a pointer of a pointer-
// to-member to avoid the ambiguity with the built-in pointer comparison operators. However,
// only pointers and pointers-to-member can be compared to Null.
*/
template< typename T >
inline bool operator!=( const T& lhs, const Null& rhs )
{
   return !rhs.equal( lhs );
}
//*************************************************************************************************

} // namespace blaze




//=================================================================================================
//
//  NULL DEFINITION
//
//=================================================================================================

#if BLAZE_USE_NULL

#undef NULL

//*************************************************************************************************
/*!\brief Global NULL pointer.
// \ingroup util
//
// This instance of the Null class replaces the NULL macro to ensure a type-safe NULL pointer.
*/
const blaze::Null NULL;
//*************************************************************************************************

#endif

#endif
