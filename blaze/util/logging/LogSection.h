//=================================================================================================
/*!
//  \file blaze/util/logging/LogSection.h
//  \brief Header file for the LogSection class
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

#ifndef _BLAZE_UTIL_LOGGING_LOGSECTION_H_
#define _BLAZE_UTIL_LOGGING_LOGSECTION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <new>
#include <sstream>
#include <blaze/system/Logging.h>
#include <blaze/util/logging/LogLevel.h>


namespace blaze {

namespace logging {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Logging section for (non-)MPI-parallel environments.
// \ingroup logging
//
// The LogSection class is an auxiliary helper class for all logging section macros. It is
// implemented as a wrapper around the Logger class and is responsible for the atomicity of
// the logging operations and for formatting any message that is written into the log file(s).
*/
class LogSection
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
          LogSection( LogLevel level );
   inline LogSection( const LogSection& ls );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~LogSection();
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   inline operator bool() const;
   inline operator std::ostream&();
   //@}
   //**********************************************************************************************

   //**Logging functions***************************************************************************
   /*!\name Logging functions */
   //@{
   template< typename Type > inline void log   ( const Type& message );
                                    void commit();
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   LogLevel          level_;    //!< The logging level of the log section.
   std::stringstream message_;  //!< Intermediate buffer for log messages.
   //@}
   //**********************************************************************************************

   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   LogSection& operator=( const LogSection& );

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
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The copy constructor for LogSection.
//
// \param ls The log section to be copied.
//
// The copy constructor is explicitly defined in order to enable its use in the log sections
// despite the non-copyable stringstream member variable.
*/
inline LogSection::LogSection( const LogSection& ls )
   : level_( ls.level_ )  // The logging level of the log section
{}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion operator to \a bool.
//
// The conversion operator returns \a true to indicate that the logging section is active.
*/
inline LogSection::operator bool() const
{
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to \a std::ostream&.
//
// The conversion operator returns a reference to a standard output stream.
*/
inline LogSection::operator std::ostream&()
{
   return message_;
}
//*************************************************************************************************




//=================================================================================================
//
//  LOGGING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Logs the given message to the log file.
//
// \param message The log message to be logged.
// \return void
*/
template< typename Type >  // Type of the log message
inline void LogSection::log( const Type& message )
{
   message_ << message;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name LogSection operators */
//@{
template< typename Type >
inline LogSection& operator<<( LogSection& logsection, const Type& message );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for the LogSection class.
// \ingroup logging
//
// \param logsection Reference to the log section.
// \param message Reference to the log message.
// \return Reference to the log section.
*/
template< typename Type >  // Type of the log message
inline LogSection& operator<<( LogSection& logsection, const Type& message )
{
   logsection.log( message );
   return logsection;
}
//*************************************************************************************************

} // namespace logging

} // namespace blaze

#endif
