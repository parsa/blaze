//=================================================================================================
/*!
//  \file blaze/util/Types.h
//  \brief Header file for basic type definitions
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

#ifndef _BLAZE_UTIL_TYPES_H_
#define _BLAZE_UTIL_TYPES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstddef>
#include <boost/cstdint.hpp>


namespace blaze {

//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::size_t
// \brief Size type of the Blaze library.
// \ingroup util
*/
using std::size_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::ptrdiff_t
// \brief Pointer difference type of the Blaze library.
// \ingroup util
*/
using std::ptrdiff_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::int8_t
// \brief 8-bit signed integer type of the Blaze library.
// \ingroup util
*/
using boost::int8_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::uint8_t
// \brief 8-bit unsigned integer type of the Blaze library.
// \ingroup util
*/
using boost::uint8_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::int16_t
// \brief 16-bit signed integer type of the Blaze library.
// \ingroup util
*/
using boost::int16_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::uint16_t
// \brief 16-bit unsigned integer type of the Blaze library.
// \ingroup util
*/
using boost::uint16_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::int32_t
// \brief 32-bit signed integer type of the Blaze library.
// \ingroup util
*/
using boost::int32_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::uint32_t
// \brief 32-bit unsigned integer type of the Blaze library.
// \ingroup util
*/
using boost::uint32_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::int64_t
// \brief 64-bit signed integer type of the Blaze library.
// \ingroup util
*/
#ifndef BOOST_NO_INT64_T
using boost::int64_t;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::uint64_t
// \brief 64-bit unsigned integer type of the Blaze library.
// \ingroup util
*/
#ifndef BOOST_NO_INT64_T
using boost::uint64_t;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The largest available signed integer data type.
// \ingroup util
*/
#ifndef BOOST_NO_INT64_T
typedef int64_t  large_t;
#else
typedef int32_t  large_t;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The largest available unsigned integer data type.
// \ingroup util
*/
#ifndef BOOST_NO_INT64_T
typedef uint64_t  ularge_t;
#else
typedef uint32_t  ularge_t;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unsigned integer data type for integral IDs.
// \ingroup util
*/
typedef ularge_t  id_t;
//*************************************************************************************************

} // namespace blaze

#endif
