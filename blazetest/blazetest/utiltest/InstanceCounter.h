//=================================================================================================
/*!
//  \file blazetest/utiltest/InstanceCounter.h
//  \brief Header file for InstanceCounter class template
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

#ifndef _BLAZETEST_UTILTEST_INSTANCECOUNTER_H_
#define _BLAZETEST_UTILTEST_INSTANCECOUNTER_H_


namespace blazetest {

namespace utiltest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for classes requiring an instance counter.
//
// The InstanceCounter class template provides the functionality to count the number of instances
// created from a particular class. The functionality can be used by either deriving publicly
// (providing public access to the instance counter) or non-publicly (hiding access to the
// instance counter) from the InstanceCounter class. The InstanceCounter class template requires
// the deriving class as template argument (CRTP pattern):

   \code
   class Resource : public InstanceCounter<Resource>
   {
      // ...
   };
   \endcode
*/
template< typename T >
class InstanceCounter
{
 protected:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline InstanceCounter();
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~InstanceCounter();
   //@}
   //**********************************************************************************************

 public:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline unsigned int getCount();
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static unsigned int counter_;  // The instance counter.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The constructor of InstanceCounter.
*/
template< typename T >
inline InstanceCounter<T>::InstanceCounter()
{
   ++counter_;
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor of InstanceCounter.
*/
template< typename T >
inline InstanceCounter<T>::~InstanceCounter()
{
   --counter_;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current count of instances.
//
// \return The current count of instances.
*/
template< typename T >
inline unsigned int InstanceCounter<T>::getCount()
{
   return counter_;
}
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename T >
unsigned int InstanceCounter<T>::counter_( 0U );

} // namespace utiltest

} // namespace blazetest

#endif
