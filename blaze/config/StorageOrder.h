//=================================================================================================
/*!
//  \file blaze/config/StorageOrder.h
//  \brief Configuration of the default storage order for all matrices of the Blaze library
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
/*!\brief The default storage order for all matrices of the Blaze library.
// \ingroup config
//
// This value specifies the default storage order for all matrices of the Blaze library.
// In case no explicit storage order is specified with the according matrix type, this
// setting is used.

   \code
   // Explicit specification of the storage order => row-major matrix
   StaticMatrix<double,3UL,3UL,rowMajor> A;

   // No explicit specification of the storage order => use of the defaultStorageOrder
   StaticMatrix<double,3UL,3UL> B;
   \endcode

// Valid settings for the defaultStorageOrder are blaze::rowMajor and blaze::columnMajor.
*/
const bool defaultStorageOrder = rowMajor;
//*************************************************************************************************

} // namespace blaze
