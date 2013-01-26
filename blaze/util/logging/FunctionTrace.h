//=================================================================================================
/*!
//  \file blaze/util/logging/FunctionTrace.h
//  \brief Header file for the FunctionTrace class
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

#ifndef _BLAZE_UTIL_LOGGING_FUNCTIONTRACE_H_
#define _BLAZE_UTIL_LOGGING_FUNCTIONTRACE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <new>
#include <string>
#include <blaze/system/Logging.h>
#include <blaze/system/Signature.h>
#include <blaze/util/NonCopyable.h>


namespace blaze {

namespace logging {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief RAII object for function tracing.
// \ingroup logging
//
// The FunctionTrace class is a 

// The FunctionTrace class is an auxiliary helper class for the tracing of function calls. It
// is implemented as a wrapper around the Logger class and is responsible for the atomicity of
// the logging operations of trace information.
*/
class FunctionTrace : private NonCopyable
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   FunctionTrace( const std::string& file, const std::string& function );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~FunctionTrace();
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string file_;      //!< The file name the traced function is contained in.
   std::string function_;  //!< The name of the traced function.
   //@}
   //**********************************************************************************************

   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   void* operator new  ( std::size_t ) throw( std::bad_alloc );
   void* operator new[]( std::size_t ) throw( std::bad_alloc );
   void* operator new  ( std::size_t, const std::nothrow_t& ) throw();
   void* operator new[]( std::size_t, const std::nothrow_t& ) throw();

   void operator delete  ( void* ) throw();
   void operator delete[]( void* ) throw();
   void operator delete  ( void*, const std::nothrow_t& ) throw();
   void operator delete[]( void*, const std::nothrow_t& ) throw();
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  BLAZE_FUNCTION_TRACE MACRO
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Function trace macro.
// \ingroup logging
//
// This macro can be used to reliably trace function calls. In case function tracing is
// activated, the traces are logged via the Logger class. The following, short example
// demonstrates how the function trace macro is used:

   \code
   int main( int argc, char** argv )
   {
      BLAZE_FUNCTION_TRACE;
      
      // ...
   }
   \endcode

// The macro should be used as the very first statement inside the function in order to
// guarantee that logging the function trace is the very first and last action of the
// function call. In case function tracing is activated, the resulting log will contain
// trace information of the following form:

   \code
   [TRACE   ][000:00:00] + Entering function 'int main()' in file 'TraceDemo.cpp'
   [TRACE   ][000:00:10] - Leaving function 'int main()' in file 'TraceDemo.cpp'
   \endcode
*/
#if BLAZE_USE_FUNCTION_TRACES
#  define BLAZE_FUNCTION_TRACE \
   blaze::logging::FunctionTrace BLAZE_FUNCTION_TRACE_OBJECT( __FILE__, BLAZE_SIGNATURE )
#else
#  define BLAZE_FUNCTION_TRACE
#endif
//*************************************************************************************************

} // namespace logging

} // namespace blaze

#endif
