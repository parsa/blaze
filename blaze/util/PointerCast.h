//=================================================================================================
/*!
//  \file blaze/util/PointerCast.h
//  \brief Cast operators for pointer types
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

#ifndef _BLAZE_UTIL_POINTERCAST_H_
#define _BLAZE_UTIL_POINTERCAST_H_


namespace blaze {

//=================================================================================================
//
//  POINTER CAST OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Pointer cast operators */
//@{
template< typename To, typename From > inline To* static_pointer_cast( From* ptr );
template< typename To, typename From > inline To* dynamic_pointer_cast( From* ptr );
template< typename To, typename From > inline To* const_pointer_cast( From* ptr);
template< typename To, typename From > inline To* reinterpret_pointer_cast( From* ptr );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Static cast for pointer types.
// \ingroup util
//
// \param ptr The pointer to be cast.
// \return The casted value.
//
// The static_pointer_cast function is used exactly as the built-in static_cast operator but
// for pointer types.

   \code
   class B { ... };
   class D : public B { ... };

   B* b = new D();                      // Base pointer to a derived class object
   D* d = static_pointer_cast<D>( b );  // Static down-cast
   \endcode
*/
template< typename To, typename From >
inline To* static_pointer_cast( From* ptr )
{
   return static_cast<To*>( ptr );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Dynamic cast for pointer types.
// \ingroup util
//
// \param ptr The pointer to be cast.
// \return The casted value.
//
// The dynamic_pointer_cast function is used exactly as the built-in dynamic_cast operator but
// for pointer types. As in case with the built-in dynamic_cast 0 is returned if the runtime
// type conversion doesn't succeed.

   \code
   class B { ... };
   class D : public B { ... };

   B* b = ...;                           // Base pointer
   D* d = dynamic_pointer_cast<D>( b );  // Dynamic down-cast
   \endcode
*/
template< typename To, typename From >
inline To* dynamic_pointer_cast( From* ptr )
{
   return dynamic_cast<To*>( ptr );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Const cast for pointer types.
// \ingroup util
//
// \param ptr The pointer to be cast.
// \return The casted value.
//
// The const_pointer_cast function is used exactly as the built-in const_cast operator but
// for pointer types.

   \code
   class A { ... };

   const A* a1;                          // Pointer to a constant A object
   A* a2 = const_pointer_cast<A>( a1 );  // Const cast to a pointer to a non-constant A object
   \endcode
*/
template< typename To, typename From >
inline To* const_pointer_cast( From* ptr)
{
   return const_cast<To*>( ptr );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reinterpret cast for pointer types.
// \ingroup util
//
// \param ptr The pointer to be cast.
// \return The casted value.
//
// The reinterpret_pointer_cast function is used exactly as the built-in reinterpret_cast
// operator but for pointer types.

   \code
   class A { ... };

   unsigned char* raw = new unsigned char[ sizeof(A)*10 ];  // Allocation of raw memory
   A* a = reinterpret_pointer_cast<A>( raw );               // Reinterpretation cast
   \endcode
*/
template< typename To, typename From >
inline To* reinterpret_pointer_cast( From* ptr )
{
   return reinterpret_cast<To*>( ptr );
}
//*************************************************************************************************




//=================================================================================================
//
//  SMART POINTER CAST OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Smart pointer cast operators */
//@{
template< typename To, template<typename> class S, typename From > inline S<To> static_pointer_cast( S<From> ptr );
template< typename To, template<typename> class S, typename From > inline S<To> dynamic_pointer_cast( S<From> ptr );
template< typename To, template<typename> class S, typename From > inline S<To> const_pointer_cast( S<From> ptr);
template< typename To, template<typename> class S, typename From > inline S<To> reinterpret_pointer_cast( S<From> ptr );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Static cast for smart pointers.
// \ingroup util
//
// \param ptr The smart pointer to be cast.
// \return The casted smart pointer.
//
// The static_pointer_cast function is used exactly as the built-in static_cast operator but
// for smart pointers.

   \code
   class B { ... };
   class D : public B { ... };

   typedef SharedPtr<B> BPtr;
   typedef SharedPtr<D> DPtr;

   BPtr b = BPtr( new D() );              // Base smart pointer to a derived class object
   DPtr d = static_pointer_cast<D>( b );  // Static down-cast
   \endcode
*/
template< typename To, template<typename> class S, typename From >
inline S<To> static_pointer_cast( S<From> ptr )
{
   return S<To>( static_cast<To*>( ptr.get() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Dynamic cast for smart pointers.
// \ingroup util
//
// \param ptr The smart pointer to be cast.
// \return The casted smart pointer.
//
// The dynamic_pointer_cast function is used exactly as the built-in dynamic_cast operator
// but for smart pointers. As in case with the built-in dynamic_cast 0 is returned if the
// runtime type conversion doesn't succeed.

   \code
   class B { ... };
   class D : public B { ... };

   typedef SharedPtr<B> BPtr;
   typedef SharedPtr<D> DPtr;

   BPtr b = ...;                           // Base smart pointer
   DPtr d = dynamic_pointer_cast<D>( b );  // Dynamic down-cast
   \endcode
*/
template< typename To, template<typename> class S, typename From >
inline S<To> dynamic_pointer_cast( S<From> ptr )
{
   return S<To>( dynamic_cast<To*>( ptr.get() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Const cast for smart pointers.
// \ingroup util
//
// \param ptr The smart pointer to be cast.
// \return The casted smart pointer.
//
// The const_pointer_cast function is used exactly as the built-in const_cast operator but
// for smart pointers.

   \code
   class A { ... };

   typedef SharedPtr<A>        APtr;
   typedef SharedPtr<const A>  ConstAPtr;

   ConstAPtr a1;                           // Smart pointer to a constant A object
   APtr a2 = const_pointer_cast<A>( a1 );  // Const cast to a smart pointer to a non-constant A object
   \endcode
*/
template< typename To, template<typename> class S, typename From >
inline S<To> const_pointer_cast( S<From> ptr )
{
   return S<To>( const_cast<To*>( ptr.get() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reinterpret cast for smart pointers.
// \ingroup util
//
// \param ptr The smart pointer to be cast.
// \return The casted smart pointer.
//
// The reinterpret_pointer_cast function is used exactly as the built-in reinterpret_cast
// operator but for smart pointers.
*/
template< typename To, template<typename> class S, typename From >
inline S<To> reinterpret_pointer_cast( S<From> ptr )
{
   return S<To>( reinterpret_cast<To*>( ptr.get() ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
