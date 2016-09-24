//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsConstructible.h
//  \brief Header file for the IsConstructible type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISCONSTRUCTIBLE_H_
#define _BLAZE_UTIL_TYPETRAITS_ISCONSTRUCTIBLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <type_traits>
#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS ISCONSTRUCTIBLE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// The IsConstructible type trait tests whether the expression

   \code
   T obj( std::declval<Args>()... );
   \endcode

// is well formed. If an object of type \a T can be created in this way, the \a value member
// constant is set to \a true, the nested type definition \a Type is set to \a TrueType, and
// the class derives from \a TrueType. Otherwise \a value is set to \a false, \a Type is
// \a FalseType and the class derives from \a FalseType.
*/
template< typename T, typename... Args >
struct IsConstructible
   : public BoolConstant< std::is_constructible<T,Args...>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the IsConstructible type trait.
// \ingroup type_traits
//
// The IsConstructible_ alias declaration provides a convenient shortcut to access the nested
// \a value of the IsConstructible class template. For instance, given the types \a T1, \a T2
// and \a T3 the following two statements are identical:

   \code
   constexpr bool value1 = IsConstructible<T1,T2,T3>::value;
   constexpr bool value2 = IsConstructible_<T1,T2,T3>;
   \endcode
*/
template< typename T, typename... Args >
constexpr bool IsConstructible_ = IsConstructible<T>::value;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS ISNOTHROWCONSTRUCTIBLE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// The IsNothrowConstructible type trait tests whether the expression

   \code
   T obj( std::declval<Args>()... );
   \endcode

// is well formed and guaranteed to not throw an exception (i.e. noexcept). If an object of type
// \a T can be created in this way, the \a value member constant is set to \a true, the nested
// type definition \a Type is set to \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to \a false, \a Type is \a FalseType and the class derives from \a FalseType.
*/
template< typename T, typename... Args >
struct IsNothrowConstructible
   : public BoolConstant< std::is_nothrow_constructible<T,Args...>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the IsNothrowConstructible type trait.
// \ingroup type_traits
//
// The IsNothrowConstructible_ alias declaration provides a convenient shortcut to access the
// nested \a value of the IsNothrowConstructible class template. For instance, given the types
// \a T1, \a T2 and \a T3 the following two statements are identical:

   \code
   constexpr bool value1 = IsNothrowConstructible<T1,T2,T3>::value;
   constexpr bool value2 = IsNothrowConstructible_<T1,T2,T3>;
   \endcode
*/
template< typename T, typename... Args >
constexpr bool IsNothrowConstructible_ = IsNothrowConstructible<T>::value;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS ISDEFAULTCONSTRUCTIBLE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// The IsDefaultConstructible type trait tests whether the expression

   \code
   T obj;
   \endcode

// is well formed. If an object of type \a T can be default constructed, the \a value member
// constant is set to \a true, the nested type definition \a Type is set to \a TrueType, and
// the class derives from \a TrueType. Otherwise \a value is set to \a false, \a Type is
// \a FalseType and the class derives from \a FalseType.
*/
template< typename T >
struct IsDefaultConstructible
   : public BoolConstant< std::is_default_constructible<T>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the IsDefaultConstructible type trait.
// \ingroup type_traits
//
// The IsDefaultConstructible_ alias declaration provides a convenient shortcut to access the
// nested \a value of the IsDefaultConstructible class template. For instance, given the type
// \a T the following two statements are identical:

   \code
   constexpr bool value1 = IsDefaultConstructible<T>::value;
   constexpr bool value2 = IsDefaultConstructible_<T>;
   \endcode
*/
template< typename T >
constexpr bool IsDefaultConstructible_ = IsDefaultConstructible<T>::value;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS ISNOTHROWDEFAULTCONSTRUCTIBLE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// The IsNothrowDefaultConstructible type trait tests whether the expression

   \code
   T obj;
   \endcode

// is well formed and guaranteed to not throw an exception (i.e. noexcept). If an object of type
// \a T can be default constructed, the \a value member constant is set to \a true, the nested
// type definition \a Type is set to \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to \a false, \a Type is \a FalseType and the class derives from \a FalseType.
*/
template< typename T >
struct IsNothrowDefaultConstructible
   : public BoolConstant< std::is_nothrow_default_constructible<T>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the IsNothrowDefaultConstructible type trait.
// \ingroup type_traits
//
// The IsNothrowDefaultConstructible_ alias declaration provides a convenient shortcut to access
// the nested \a value of the IsNothrowDefaultConstructible class template. For instance, given
// the type \a T the following two statements are identical:

   \code
   constexpr bool value1 = IsNothrowDefaultConstructible<T>::value;
   constexpr bool value2 = IsNothrowDefaultConstructible_<T>;
   \endcode
*/
template< typename T >
constexpr bool IsNothrowDefaultConstructible_ = IsNothrowDefaultConstructible<T>::value;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS ISCOPYCONSTRUCTIBLE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// The IsCopyConstructible type trait tests whether the expression

   \code
   T obj( std::declval<T>() );
   \endcode

// is well formed. If an object of type \a T can be copy constructed, the \a value member
// constant is set to \a true, the nested type definition \a Type is set to \a TrueType, and
// the class derives from \a TrueType. Otherwise \a value is set to \a false, \a Type is
// \a FalseType and the class derives from \a FalseType.
*/
template< typename T >
struct IsCopyConstructible
   : public BoolConstant< std::is_copy_constructible<T>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the IsCopyConstructible type trait.
// \ingroup type_traits
//
// The IsCopyConstructible_ alias declaration provides a convenient shortcut to access the nested
// \a value of the IsCopyConstructible class template. For instance, given the type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = IsCopyConstructible<T>::value;
   constexpr bool value2 = IsCopyConstructible_<T>;
   \endcode
*/
template< typename T >
constexpr bool IsCopyConstructible_ = IsCopyConstructible<T>::value;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS ISNOTHROWCOPYCONSTRUCTIBLE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// The IsNothrowCopyConstructible type trait tests whether the expression

   \code
   T obj( std::declval<T>() );
   \endcode

// is well formed and guaranteed to not throw an exception (i.e. noexcept). If an object of type
// \a T can be copy constructed, the \a value member constant is set to \a true, the nested type
// definition \a Type is set to \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to \a false, \a Type is \a FalseType and the class derives from \a FalseType.
*/
template< typename T >
struct IsNothrowCopyConstructible
   : public BoolConstant< std::is_copy_constructible<T>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the IsNothrowCopyConstructible type trait.
// \ingroup type_traits
//
// The IsNothrowCopyConstructible_ alias declaration provides a convenient shortcut to access the
// nested \a value of the IsNothrowCopyConstructible class template. For instance, given the type
// \a T the following two statements are identical:

   \code
   constexpr bool value1 = IsNothrowCopyConstructible<T>::value;
   constexpr bool value2 = IsNothrowCopyConstructible_<T>;
   \endcode
*/
template< typename T >
constexpr bool IsNothrowCopyConstructible_ = IsNothrowCopyConstructible<T>::value;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS ISMOVECONSTRUCTIBLE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// The IsMoveConstructible type trait tests whether the expression

   \code
   T obj( std::move( std::declval<T>() ) );
   \endcode

// is well formed. If an object of type \a T can be move constructed, the \a value member
// constant is set to \a true, the nested type definition \a Type is set to \a TrueType, and
// the class derives from \a TrueType. Otherwise \a value is set to \a false, \a Type is
// \a FalseType and the class derives from \a FalseType.
*/
template< typename T >
struct IsMoveConstructible
   : public BoolConstant< std::is_move_constructible<T>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the IsMoveConstructible type trait.
// \ingroup type_traits
//
// The IsMoveConstructible_ alias declaration provides a convenient shortcut to access the nested
// \a value of the IsMoveConstructible class template. For instance, given the type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = IsMoveConstructible<T>::value;
   constexpr bool value2 = IsMoveConstructible_<T>;
   \endcode
*/
template< typename T >
constexpr bool IsMoveConstructible_ = IsMoveConstructible<T>::value;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS ISNOTHROWMOVECONSTRUCTIBLE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// The IsNothrowMoveConstructible type trait tests whether the expression

   \code
   T obj( std::move( std::declval<T>() ) );
   \endcode

// is well formed and guaranteed to not throw an exception (i.e. noexcept). If an object of type
// \a T can be move constructed, the \a value member constant is set to \a true, the nested type
// definition \a Type is set to \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to \a false, \a Type is \a FalseType and the class derives from \a FalseType.
*/
template< typename T >
struct IsNothrowMoveConstructible
   : public BoolConstant< std::is_nothrow_move_constructible<T>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the IsNothrowMoveConstructible type trait.
// \ingroup type_traits
//
// The IsNothrowMoveConstructible_ alias declaration provides a convenient shortcut to access the
// nested \a value of the IsNothrowMoveConstructible class template. For instance, given the type
// \a T the following two statements are identical:

   \code
   constexpr bool value1 = IsNothrowMoveConstructible<T>::value;
   constexpr bool value2 = IsNothrowMoveConstructible_<T>;
   \endcode
*/
template< typename T >
constexpr bool IsNothrowMoveConstructible_ = IsNothrowMoveConstructible<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
