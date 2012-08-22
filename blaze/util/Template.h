//=================================================================================================
/*!
//  \file blaze/util/Template.h
//  \brief Header file for nested template disabiguation
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

#ifndef _BLAZE_UTIL_TEMPLATE_H_
#define _BLAZE_UTIL_TEMPLATE_H_


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compiler specific patch for nested template disambiguation.
// \ingroup util
//
// The BLAZE_TEMPLATE is a patch for the Microsoft Visual C++ compiler that does not correctly
// parse definitions of nested templates of the following form:

   \code
   template< typename T >
   class Alloc {
    public:
      ...
      template< typename Other >
      class rebind {
       public:
         typedef Alloc<Other> other;
      };
      ...
   };

   typedef Alloc<int>  AI;
   typedef AI::template rebind<double>::other  Other;  // Compilation error with Visual C++
   \endcode

// In order to circumvent this compilation error, the BLAZE_TEMPLATE macro should be used
// instead the \a template keyword:

   \code
   ...
   typedef AI::BLAZE_TEMPLATE rebind<double>::other  Other;  // No compilation errors
   \endcode
*/
#if defined(_MSC_VER)
#  define BLAZE_TEMPLATE
#else
#  define BLAZE_TEMPLATE template
#endif
/*! \endcond */
//*************************************************************************************************

#endif
