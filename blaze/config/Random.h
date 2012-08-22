//=================================================================================================
/*!
//  \file blaze/config/Random.h
//  \brief Configuration of the random number generator used in the Blaze library
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
/*!\brief Type of the random number generator of the Blaze library.
// \ingroup config
//
// This type definition represents the type of the random number generated used in the Blaze
// library. The default random number generator is the boost::mt19937 mersenne-twister pseudo
// random number generator. For more information see the class description of the boost library:
//
//   http://www.boost.org/doc/libs/1_35_0/libs/random/random-generators.html#mersenne_twister\n
//   http://www.boost.org/doc/libs/1_35_0/boost/random/mersenne_twister.hpp
*/
typedef boost::mt19937  RNG;
//*************************************************************************************************

} // namespace blaze
