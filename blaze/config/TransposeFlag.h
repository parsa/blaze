//=================================================================================================
/*!
//  \file blaze/config/TransposeFlag.h
//  \brief Configuration of the default transpose flag for all vectors of the Blaze library
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


namespace blaze {

//*************************************************************************************************
/*!\brief The default transpose flag for all vectors of the Blaze library.
// \ingroup config
//
// This value specifies the default transpose flag for all vector of the Blaze library.
// In case no explicit transpose flag is specified with the according vector type, this
// setting is used.

   \code
   // Explicit specification of the transpose flag => column vector
   StaticVector<double,3UL,columnMajor> a;

   // No explicit specification of the transpose flag => use of the defaultTransposeFlag
   StaticVector<double,3UL> b;
   \endcode

// Valid settings for the defaultTransposeFlag are blaze::rowVector and blaze::columnVector.
*/
const bool defaultTransposeFlag = columnVector;
//*************************************************************************************************

} // namespace blaze
