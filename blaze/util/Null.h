//=================================================================================================
/*!
//  \file blaze/util/Null.h
//  \brief Header file for a safe C++ NULL pointer implementation
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
*/
//=================================================================================================

#ifndef _BLAZE_UTIL_NULL_H_
#define _BLAZE_UTIL_NULL_H_


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

#ifdef NULL
#  undef NULL
#endif

//*************************************************************************************************
/*!\brief Global NULL pointer.
// \ingroup util
//
// This instance of the Null class replaces the NULL macro to ensure a type-safe NULL pointer.
*/
const blaze::Null NULL;
//*************************************************************************************************

#endif
