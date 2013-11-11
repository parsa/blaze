//=================================================================================================
/*!
//  \file blazemark/util/Indices.h
//  \brief Header file for the Indices class
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
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
