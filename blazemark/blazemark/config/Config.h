//=================================================================================================
/*!
//  \file blazemark/config/Config.h
//  \brief General configuration file for the blaze benchmark suite
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


//=================================================================================================
//
//  BENCHMARK CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
const bool runClike     ( true );  //!< C-like, manual implementation
const bool runClassic   ( true );  //!< Classical operator overloading benchmark tests
const bool runBLAS      ( true );  //!< BLAS benchmark tests
const bool runBlaze     ( true );  //!< Blaze benchmark tests
const bool runBoost     ( true );  //!< Boost uBLAS benchmark tests
const bool runBlitz     ( true );  //!< Blitz++ benchmark tests
const bool runGMM       ( true );  //!< GMM++ benchmark tests
const bool runArmadillo ( true );  //!< Armadillo benchmark tests
const bool runFLENS     ( true );  //!< FLENS benchmark tests
const bool runMTL       ( true );  //!< MTL benchmark tests
const bool runEigen     ( true );  //!< Eigen benchmark tests
//*************************************************************************************************


//*************************************************************************************************
const size_t reps     ( 3     );  //!< Configuration of the number of benchmark repetitions
const double runtime  ( 2.0   );  //!< Target runtime for a benchmark measurement
const double maxtime  ( 600.0 );  //!< Maximum runtime of a single benchmark measurement [s]
//*************************************************************************************************


//*************************************************************************************************
const double deviation( 5.0 );  //!< Maximum allowed deviation of the average benchmark time from the minimum time
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Element type for all benchmarks.
//
// This type definition specifies the data type used as element type for all dense and sparse
// vectors and matrices in all benchmarks. Possible settings for this type are <a>char, signed
// char, unsigned char, wchar_t, short, unsigned short, int, unsigned int, long, unsigned long,
// float, double, long double and complex</a> with the previously named types as element type.
// The default setting is \a double.
*/
typedef double  element_t;
//*************************************************************************************************





//=================================================================================================
//
//  RANDOM NUMBER CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
const size_t seed( 128753984 );  //!< Seed for the random number generator.
//*************************************************************************************************




//=================================================================================================
//
//  EIGEN CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Index data type for the Eigen sparse matrices.
//
// This type definition offers the possibility to specify the index type for the sparse matrices
// in all Eigen benchmarks.
*/
typedef int  EigenSparseIndexType;
//*************************************************************************************************
