//=================================================================================================
/*!
//  \file blazemark/util/Benchmarks.h
//  \brief Header file for the Benchmarks class
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

#ifndef _BLAZEMARK_UTIL_BENCHMARKS_H_
#define _BLAZEMARK_UTIL_BENCHMARKS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstring>
#include <sstream>
#include <stdexcept>
#include <blazemark/system/Config.h>


namespace blazemark {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Configuration data structure for the selection of benchmarks.
//
// This auxiliary data structure represents the selection of different benchmarks (Blaze, Boost,
// Blitz, ...) for a benchmark run.
*/
struct Benchmarks
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Benchmarks();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Copy assignment operator********************************************************************
   // No explicitly declared copy assignment operator.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   bool runClike;      //!< Flag value for the C-like benchmark kernels.
                       /*!< In case the runClike flag is set to \a true and in case a C-like kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runClike flag is set to \a false, the C-like
                            kernel will be skipped.*/
   bool runClassic;    //!< Flag value for the classic benchmark kernels.
                       /*!< In case the runClassic flag is set to \a true and in case a classic kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runClassic flag is set to \a false, the
                            classic kernel will be skipped.*/
   bool runBLAS;       //!< Flag value of the BLAS benchmark kernels.
                       /*!< In case the runBLAS flag is set to \a true and in case a BLAS kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runBLAS flag is set to \a false, the BLAS
                            kernel will be skipped.*/
   bool runBlaze;      //!< Flag value for the Blaze benchmark kernels.
                       /*!< In case the runBlaze flag is set to \a true and in case a Blaze kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runBlaze flag is set to \a false, the
                            Blaze kernel will be skipped.*/
   bool runBoost;      //!< Flag value for the Boost uBLAS benchmark kernels.
                       /*!< In case the runBlaze flag is set to \a true and in case a Boost kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runBoost flag is set to \a false, the
                            Boost kernel will be skipped.*/
   bool runBlitz;      //!< Flag value for the Blitz benchmark kernels.
                       /*!< In case the runBlaze flag is set to \a true and in case a Blitz kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runBlitz flag is set to \a false, the
                            Blitz kernel will be skipped.*/
   bool runGMM;        //!< Flag value for the GMM benchmark kernels.
                       /*!< In case the runGMM flag is set to \a true and in case a GMM kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runGMM flag is set to \a false, the
                            GMM kernel will be skipped.*/
   bool runArmadillo;  //!< Flag value for the Armadillo benchmark kernels.
                       /*!< In case the runArmadillo flag is set to \a true and in case an Armadillo
                            kernel is available for a particular benchmark, the kernel is included
                            in the benchmark tests. In case the runArmadillo flag is set to \a false,
                            the Armadillo kernel will be skipped.*/
   bool runFLENS;      //!< Flag value for the FLENS benchmark kernels.
                       /*!< In case the runFLENS flag is set to \a true and in case a FLENS kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runFlENS flag is set to \a false, the
                            FLENS kernel will be skipped. */
   bool runMTL;        //!< Flag value for the MTL benchmark kernels.
                       /*!< In case the runMTL flag is set to \a true and in case a MTL kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runMTL flag is set to \a false, the
                            MTL kernel will be skipped.*/
   bool runEigen;      //!< Flag value for the Eigen benchmark kernels.
                       /*!< In case the runBlaze flag is set to \a true and in case a Eigen kernel
                            is available for a particular benchmark, the kernel is included in the
                            benchmark tests. In case the runEigen flag is set to \a false, the
                            Eigen kernel will be skipped.*/
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default constructor for the Benchmarks class.
*/
inline Benchmarks::Benchmarks()
   : runClike    ( blazemark::runClike     )  // Flag value for the C-like benchmark kernels
   , runClassic  ( blazemark::runClassic   )  // Flag value for the classic benchmark kernels
   , runBLAS     ( blazemark::runBLAS      )  // Flag value for the BLAS benchmark kernels
   , runBlaze    ( blazemark::runBlaze     )  // Flag value for the Blaze benchmark kernels
   , runBoost    ( blazemark::runBoost     )  // Flag value for the Boost uBLAS benchmark kernels
   , runBlitz    ( blazemark::runBlitz     )  // Flag value for the Blitz benchmark kernels
   , runGMM      ( blazemark::runGMM       )  // Flag value for the GMM benchmark kernels
   , runArmadillo( blazemark::runArmadillo )  // Flag value for the Armadillo benchmark kernels
   , runFLENS    ( blazemark::runFLENS     )  // Flag value for the FLENS benchmark kernels
   , runMTL      ( blazemark::runMTL       )  // Flag value for the MTL benchmark kernels
   , runEigen    ( blazemark::runEigen     )  // Flag value for the Eigen benchmark kernels
{}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Parsing the command line arguments to configure the given benchmarks data structure.
//
// \param argc The total number of command line arguments.
// \param argv The array of command line arguments.
// \param benchmarks The benchmark data structure to be configured.
// \return void
// \exception std::invalid_argument Unknown command line argument.
//
// This function parses the command line arguments to configure the given benchmarks data
// structure. The following command line options will be recognized:
//
//   - \a -clike: Activates the C-like kernels.
//   - \a -no-clike: Deactivates the C-like kernels.
//   - \a -only-clike: Activates the C-like kernels and deactivates all other.
//   - \a -classic: Activates the classic kernels.
//   - \a -no-classic: Deactivates the classic kernels.
//   - \a -only-classic: Activates the classic kernels and deactivates all other.
//   - \a -blas: Activates the BLAS kernels.
//   - \a -no-blas: Deactivates the BLAS kernels.
//   - \a -only-blas: Activates the BLAS kernels and deactivates all other.
//   - \a -blaze: Activates the Blaze kernels.
//   - \a -no-blaze: Deactivates the Blaze kernels.
//   - \a -only-blaze: Activates the Blaze kernels and deactivates all other.
//   - \a -boost: Activates the Boost kernels.
//   - \a -no-boost: Deactivates the Boost kernels.
//   - \a -only-boost: Activates the Boost kernels and deactivates all other.
//   - \a -blitz: Activates the Blitz kernels.
//   - \a -no-blitz: Deactivates the Blitz kernels.
//   - \a -only-blitz: Activates the Blitz kernels and deactivates all other.
//   - \a -gmm: Activates the GMM kernels.
//   - \a -no-gmm: Deactivates the GMM kernels.
//   - \a -only-gmm: Activates the GMM kernels and deactivates all other.
//   - \a -armadillo: Activates the Armadillo kernels.
//   - \a -no-armadillo: Deactivates the Armadillo kernels.
//   - \a -only-armadillo: Activates the Armadillo kernels and deactivates all other.
//   - \a -flens: Activates the FLENS kernels.
//   - \a -no-flens: Deactivates the FLENS kernels.
//   - \a -only-flens: Activates the FLENS kernels and deactivates all other.
//   - \a -mtl: Activates the MTL kernels.
//   - \a -no-mtl: Deactivates the MTL kernels.
//   - \a -only-mtl: Activates the MTL kernels and deactivates all other.
//   - \a -eigen: Activates the Eigen kernels.
//   - \a -no-eigen: Deactivates the Eigen kernels.
//   - \a -only-eigen: Activates the Eigen kernels and deactivates all other.
//
// In case an unknown command line option is encountered, a \a std::invalid_argument exception
// is thrown.
*/
inline void parseCommandLineArguments( int argc, char** argv, Benchmarks& benchmarks )
{
   for( int i=1; i<argc; ++i )
   {
      if( std::strcmp( argv[i], "-clike" ) == 0 ) {
         benchmarks.runClike = true;
      }
      else if( std::strcmp( argv[i], "-no-clike" ) == 0 ) {
         benchmarks.runClike = false;
      }
      else if( std::strcmp( argv[i], "-only-clike" ) == 0 ) {
         benchmarks.runClike     = true;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-classic" ) == 0 ) {
         benchmarks.runClassic = true;
      }
      else if( std::strcmp( argv[i], "-no-classic" ) == 0 ) {
         benchmarks.runClassic = false;
      }
      else if( std::strcmp( argv[i], "-only-classic" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = true;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-blas" ) == 0 ){
         benchmarks.runBLAS = true;
      }
      else if( std::strcmp( argv[i], "-no-blas" ) == 0 ) {
         benchmarks.runBLAS = false;
      }
      else if( std::strcmp( argv[i], "-only-blas" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = true;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-blaze" ) == 0 ) {
         benchmarks.runBlaze = true;
      }
      else if( std::strcmp( argv[i], "-no-blaze" ) == 0 ) {
         benchmarks.runBlaze = false;
      }
      else if( std::strcmp( argv[i], "-only-blaze" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = true;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-boost" ) == 0 ) {
         benchmarks.runBoost = true;
      }
      else if( std::strcmp( argv[i], "-no-boost" ) == 0 ) {
         benchmarks.runBoost = false;
      }
      else if( std::strcmp( argv[i], "-only-boost" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = true;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-blitz" ) == 0 ) {
         benchmarks.runBlitz = true;
      }
      else if( std::strcmp( argv[i], "-no-blitz" ) == 0 ) {
         benchmarks.runBlitz = false;
      }
      else if( std::strcmp( argv[i], "-only-blitz" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = true;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-gmm" ) == 0 ) {
         benchmarks.runGMM = true;
      }
      else if( std::strcmp( argv[i], "-no-gmm" ) == 0 ) {
         benchmarks.runGMM = false;
      }
      else if( std::strcmp( argv[i], "-only-gmm" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = true;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-armadillo" ) == 0 ) {
         benchmarks.runArmadillo = true;
      }
      else if( std::strcmp( argv[i], "-no-armadillo" ) == 0 ) {
         benchmarks.runArmadillo = false;
      }
      else if( std::strcmp( argv[i], "-only-armadillo" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = true;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-flens" ) == 0 ) {
         benchmarks.runFLENS = true;
      }
      else if( std::strcmp( argv[i], "-no-flens" ) == 0 ) {
         benchmarks.runFLENS = false;
      }
      else if( std::strcmp( argv[i], "-only-flens" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = true;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-mtl" ) == 0 ) {
         benchmarks.runMTL = true;
      }
      else if( std::strcmp( argv[i], "-no-mtl" ) == 0 ) {
         benchmarks.runMTL = false;
      }
      else if( std::strcmp( argv[i], "-only-mtl" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = true;
         benchmarks.runEigen     = false;
      }
      else if( std::strcmp( argv[i], "-eigen" ) == 0 ) {
         benchmarks.runEigen = true;
      }
      else if( std::strcmp( argv[i], "-no-eigen" ) == 0 ) {
         benchmarks.runEigen = false;
      }
      else if( std::strcmp( argv[i], "-only-eigen" ) == 0 ) {
         benchmarks.runClike     = false;
         benchmarks.runClassic   = false;
         benchmarks.runBLAS      = false;
         benchmarks.runBlaze     = false;
         benchmarks.runBoost     = false;
         benchmarks.runBlitz     = false;
         benchmarks.runGMM       = false;
         benchmarks.runArmadillo = false;
         benchmarks.runFLENS     = false;
         benchmarks.runMTL       = false;
         benchmarks.runEigen     = true;
      }
      else {
         std::ostringstream oss;
         oss << " Unknown command line argument: '" << argv[i] << "'";
         throw std::invalid_argument( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace blazemark

#endif
