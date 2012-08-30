//=================================================================================================
/*!
//  \file src/blaze/DVecTSVecMult.cpp
//  \brief Source file for the dense vector/sparse vector outer product benchmark
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Functions.h>
#include <blaze/math/Infinity.h>
#include <blaze/util/Timing.h>
#include <blazemark/blaze/DVecTSVecMult.h>
#include <blazemark/boost/DVecTSVecMult.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Benchmarks.h>
#include <blazemark/util/Indices.h>
#include <blazemark/util/Parser.h>
#include <blazemark/util/SparseRun.h>


//*************************************************************************************************
// Using declarations
//*************************************************************************************************

using blazemark::Benchmarks;
using blazemark::Parser;
using blazemark::SparseRun;




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Estimating the necessary number of steps for each benchmark.
//
// \param run The parameters for the benchmark run.
// \return void
//
// This function estimates the necessary number of steps for the given benchmark based on the
// performance of the Blaze library.
*/
void estimateSteps( SparseRun& run )
{
   using blazemark::element_t;
   using blaze::rowVector;
   using blaze::columnVector;
   using blaze::rowMajor;

   const size_t N( run.getSize() );
   const size_t F( run.getNonZeros() );

   blaze::DynamicVector<element_t,columnVector> a( N, 0.1 );
   blaze::CompressedVector<element_t,rowVector> b( N, F );
   blaze::CompressedMatrix<element_t,rowMajor> A( N, N );
   blaze::timing::WcTimer timer;
   double wct( 0.0 );
   size_t steps( 1UL );

   {
      blazemark::Indices indices( N, F );
      for( blazemark::Indices::Iterator it=indices.begin(); it!=indices.end(); ++it ) {
         b[*it] = element_t(0.1);
      }
   }

   while( true ) {
      timer.start();
      for( size_t i=0UL; i<steps; ++i ) {
         A = a * b;
      }
      timer.end();
      wct = timer.last();
      if( wct >= 0.2 ) break;
      steps *= 2UL;
   }

   if( A.rows() != N )
      std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

   run.setSteps( blaze::max( 1UL, ( blazemark::runtime * steps ) / timer.last() ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  BENCHMARK FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Dense vector/sparse vector outer product benchmark function.
//
// \param runs The specified benchmark runs.
// \param benchmarks The selection of benchmarks.
// \return void
*/
void dvectsvecmult( std::vector<SparseRun>& runs, Benchmarks benchmarks )
{
   std::cout << std::left;

   std::sort( runs.begin(), runs.end() );

   size_t slowSize( blaze::inf );
   for( std::vector<SparseRun>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
      if( run->getSteps() == 0UL ) {
         if( run->getSize() < slowSize ) {
            estimateSteps( *run );
            if( run->getSteps() == 1UL )
               slowSize = run->getSize();
         }
         else run->setSteps( 1UL );
      }
   }

   if( benchmarks.runBlaze ) {
      std::vector<SparseRun>::iterator run=runs.begin();
      while( run != runs.end() ) {
         const float fill( run->getFillingDegree() );
         std::cout << "   Blaze (" << fill << "% filled) [MFlop/s]:\n";
         for( ; run!=runs.end(); ++run ) {
            if( run->getFillingDegree() != fill ) break;
            const size_t N    ( run->getSize()     );
            const size_t F    ( run->getNonZeros() );
            const size_t steps( run->getSteps()    );
            run->setBlazeResult( blazemark::blaze::dvectsvecmult( N, F, steps ) );
            const double mflops( ( N*F ) * steps / run->getBlazeResult() / 1E6 );
            std::cout << "     " << std::setw(12) << N << mflops << std::endl;
         }
      }
   }

   if( benchmarks.runBoost ) {
      std::vector<SparseRun>::iterator run=runs.begin();
      while( run != runs.end() ) {
         const float fill( run->getFillingDegree() );
         std::cout << "   Boost uBLAS (" << fill << "% filled) [MFlop/s]:\n";
         for( ; run!=runs.end(); ++run ) {
            if( run->getFillingDegree() != fill ) break;
            const size_t N    ( run->getSize()     );
            const size_t F    ( run->getNonZeros() );
            const size_t steps( run->getSteps()    );
            run->setBoostResult( blazemark::boost::dvectsvecmult( N, F, steps ) );
            const double mflops( ( N*F ) * steps / run->getBoostResult() / 1E6 );
            std::cout << "     " << std::setw(12) << N << mflops << std::endl;
         }
      }
   }

   for( std::vector<SparseRun>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
      std::cout << *run;
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The main function for the dense vector/sparse vector outer product benchmark.
//
// \param argc The total number of command line arguments.
// \param argv The array of command line arguments.
// \return void
*/
int main( int argc, char** argv )
{
   std::cout << "\n Dense Vector/Sparse Vector Outer Product:\n";

   Benchmarks benchmarks;
   parseCommandLineArguments( argc, argv, benchmarks );

   const std::string installPath( INSTALL_PATH );
   const std::string parameterFile( installPath + "/params/dvectsvecmult.prm" );
   Parser<SparseRun> parser;
   std::vector<SparseRun> runs;

   try {
      parser.parse( parameterFile.c_str(), runs );
   }
   catch( std::exception& ex ) {
      std::cerr << "   Error during parameter extraction: " << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   try {
      dvectsvecmult( runs, benchmarks );
   }
   catch( std::exception& ex ) {
      std::cerr << "   Error during benchmark execution: " << ex.what() << "\n";
      return EXIT_FAILURE;
   }
}
//*************************************************************************************************
