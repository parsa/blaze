//=================================================================================================
/*!
//  \file src/main/MemorySweep.cpp
//  \brief Source file for the Blaze memory sweep
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

#include <cstddef>
#include <cstdlib>
#include <iostream>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The main function for the Blaze memory sweep.
//
// \param argc Number of command line arguments.
// \param argv Array of command line arguments.
// \return Success code for the execution.
*/
int main( int argc, char** argv )
{
   if( argc != 2 ) {
      std::cerr << " Invalid use of program 'MemorySweep'!\n"
                << "   Use: ./memorysweep <number_of_megabytes>\n" << std::endl;
      return EXIT_FAILURE;
   }

   const std::size_t M( static_cast<std::size_t>( atoi( argv[1] ) ) );
   const std::size_t N( M * 1000000UL / sizeof(double) );

   double* v = new double[N];

   std::cout << "\n Freeing " << M << " MByte of main memory..." << std::endl;

   for( std::size_t i=0UL; i<N; ++i ) {
      if( ( i % (N/10UL) ) == 0UL )
         std::cout << "\r   Initializing the memory: " << (100.0*i)/N << "%  " << std::flush;
      v[i] = 0.0;
   }

   std::cout << "\r   Initializing the memory: 100%\n" << std::endl;

   delete [] v;
}
//*************************************************************************************************
