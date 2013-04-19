//=================================================================================================
/*!
//  \file blazetest/system/Types.h
//  \brief Type settings for the blaze test suite
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

#ifndef _BLAZETEST_SYSTEM_TYPES_H_
#define _BLAZETEST_SYSTEM_TYPES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <complex>
#include <cstddef>


namespace blazetest {

//*************************************************************************************************
/*!\class blazetest::size_t
// \brief Size type of the Blaze test suite.
*/
using std::size_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blazetest::ptrdiff_t
// \brief Pointer difference type of the Blaze test suite.
*/
using std::ptrdiff_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blazetest::complex
// \brief Complex data type of the Blaze test suite.
*/
using std::complex;
//*************************************************************************************************

} // namespace blazetest

#endif
