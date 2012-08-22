//=================================================================================================
/*!
//  \file blaze/util/threadpool/FuncWrapper.h
//  \brief Wrapper class for scheduled functions/functors
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

#ifndef _BLAZE_UTIL_THREADPOOL_FUNCWRAPPER_H_
#define _BLAZE_UTIL_THREADPOOL_FUNCWRAPPER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits.hpp>
#include <boost/utility/result_of.hpp>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/threadpool/Task.h>


namespace blaze {

namespace threadpool {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Wrapper class for scheduled functions/functors.
// \ingroup threads
//
// The FuncWrapper class is a wrapper for any callable function or functor of the following
// structure:

   \code
   // Definition of a function with no arguments that returns void
   void function() { ... }

   // Definition of a functor taking no arguments and returning void
   struct Functor
   {
      void operator()() { ... }
   };
   \endcode
*/
template< typename Callable >  // Type of the function/functor
class FuncWrapper : public Task
{
 public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline FuncWrapper( Callable func );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~FuncWrapper();
   //@}
   //**********************************************************************************************

   //**Execution functions*************************************************************************
   /*!\name Execution functions */
   //@{
   virtual void run();
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Callable func_;  //!< The wrapped function/functor.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_STATIC_ASSERT( boost::function_traits<Callable()>::arity == 0 );
   BLAZE_STATIC_ASSERT( boost::is_void< typename boost::result_of<Callable()>::type >::value );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the FuncWrapper class.
//
// \param func The given user task.
*/
template< typename Callable >  // Type of the function/functor
inline FuncWrapper<Callable>::FuncWrapper( Callable func )
   : func_( func )  // The wrapped function/functor
{}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the FuncWrapper class.
*/
template< typename Callable >  // Type of the function/functor
FuncWrapper<Callable>::~FuncWrapper()
{}
//*************************************************************************************************




//=================================================================================================
//
//  EXECUTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Executes the wrapped function/functor.
//
// \return void
*/
template< typename Callable >  // Type of the function/functor
void FuncWrapper<Callable>::run()
{
   func_();
}
//*************************************************************************************************

} // namespace threadpool

} // namespace blaze

#endif
