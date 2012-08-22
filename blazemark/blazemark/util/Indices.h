//=================================================================================================
/*!
//  \file blazemark/util/Indices.h
//  \brief Header file for the Indices class
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

#ifndef _BLAZEMARK_UTIL_INDICES_H_
#define _BLAZEMARK_UTIL_INDICES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <set>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for the generation of random indices.
//
// This auxiliary class generates random indices for sparse data structures (vector,
// matrices, etc).
*/
class Indices
{
 public:
   //**Type definitions****************************************************************************
   typedef std::set<size_t>::const_iterator  Iterator;  //!< Iterator over the generated indices.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline Indices( size_t N, size_t F );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Iterator begin() const;
   inline Iterator end  () const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::set<size_t> indices_;  //!< The generated indices for the sparse data structure.
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
/*!\brief The constructor for the Indices class.
//
// \param N The size of the sparse vector or matrix.
// \param F The filling degree of the sparse vector or matrix.
// \exception std::invalid_argument Invalid number of sparse elements.
//
// This constructor initializes an Indices object by generating \a F indices for a sparse
// data structur of size \a N. In case \a F is larger than \a N, a \a std::invalid_argument
// exception is thrown.
*/
inline Indices::Indices( size_t N, size_t F )
   : indices_()  // The generated indices for the sparse data structure
{
   if( F > N )
      throw std::invalid_argument( "Invalid number of sparse elements" );

   while( indices_.size() < F ) {
      indices_.insert( ::blaze::rand<size_t>(0,N-1) );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns an iterator to the beginning of the vector.
//
// \return Iterator to the beginning of the vector.
*/
inline Indices::Iterator Indices::begin() const
{
   return indices_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the vector.
//
// \return Iterator just past the last element of the vector.
*/
inline Indices::Iterator Indices::end() const
{
   return indices_.end();
}
//*************************************************************************************************

} // namespace blazemark

#endif
