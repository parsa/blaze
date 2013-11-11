//*************************************************************************************************
/*!
//  \file src/math/solvers/Lemke.cpp
//  \brief The Lemke pivoting algorithm for solving LCPs
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
//*************************************************************************************************


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <ostream>
#include <boost/format.hpp>
#include <blaze/math/Accuracy.h>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/Infinity.h>
#include <blaze/math/solvers/Lemke.h>
#include <blaze/util/Assert.h>
#include <blaze/util/ColorMacros.h>
#include <blaze/util/logging/DebugSection.h>
#include <blaze/util/Random.h>


namespace blaze {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default constructor for the Lemke class.
*/
Lemke::Lemke()
   : Solver()      // Initialization of the base class
   , basics_()     // TODO
   , nonbasics_()  // TODO
   , M_()          // TODO
   , Q_()          // TODO
{}
//*************************************************************************************************




//=================================================================================================
//
//  SOLVER FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief TODO
//
// \param lcp TODO
// \return void
//
// TODO
*/
bool Lemke::solve( LCP& lcp )
{
   bool converged( false );
   VecN coverVector( lcp.size(), 1 );

   size_t it( 0 );
   for( ; !converged && it<maxIterations_; ++it )
   {
      bool solved( solve( lcp, coverVector ) );

      lastPrecision_ = lcp.residual();

      if( lastPrecision_ < threshold_ )
      {
         BLAZE_LOG_DEBUG_SECTION( log ) {
            if( !solved )
               log << BLAZE_YELLOW << "      WARNING: The LCP was only solved approximately." << BLAZE_OLDCOLOR;
         }

         converged = true;
      }
      else
      {
         BLAZE_LOG_DEBUG_SECTION( log ) {
            if( solved )
               log << BLAZE_YELLOW << "      WARNING: The LCP was supposedly solved but not accurately enough (" << lastPrecision_ << ")." << BLAZE_OLDCOLOR;
            else
               log << BLAZE_YELLOW << "      WARNING: The LCP was NOT solved." << BLAZE_OLDCOLOR;
         }

         // Retrying with a different randomly chosen cover vector
         for( size_t i=0; i<coverVector.size(); ++i )
            coverVector[i] = rand<real>( 0.1, 9.91 );
      }
   }

   BLAZE_LOG_DEBUG_SECTION( log ) {
      if( converged && it == 1 )
         log << "      Solved the LCP on first try.";
      else if( converged && it > 1 )
         log << BLAZE_YELLOW << "      WARNING: Solved the LCP in " << it << " tries." << BLAZE_OLDCOLOR;
      else
         log << BLAZE_YELLOW << "      WARNING: Did not solve the LCP within accuracy. (" << lastPrecision_ << ")" << BLAZE_OLDCOLOR;
   }

   lastIterations_ = it;

   return converged;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
//
// \param lcp TODO
// \param d TODO
// \return TODO
//
// TODO
*/
bool Lemke::solve( LCP& lcp, const VecN& d )
{
   const size_t n( lcp.size() );

   const CMatMxN& A( lcp.A_ );
   const VecN&    b( lcp.b_ );

   basics_.resize( n, false );
   nonbasics_.resize( n+1, false );
   M_.resize( n, n+1, false );
   Q_.resize( n, n+1, false );

   // Merging q into Q' = [q; Q]
   for( size_t i=0; i<n; ++i )
   {
      Q_(i,0) = b[i];

      // Constructing a lexicographically positive identity matrix Q
      for( size_t j=1; j<Q_.columns(); ++j ) {
         if( j == i+1 )
            Q_(i,j) = 1;
         else
            Q_(i,j) = 0;
      }
   }

   // Preparing the augmented LCP with M' = [d; M]
   M_ = real( 0 );
   for( size_t i=0; i<n; ++i ) {
      M_(i,0) = d[i];

      for( CMatMxN::ConstIterator element=A.begin(i); element!=A.end(i); ++element ) {
         M_(i,element->index()+1) = element->value();
      }
   }

   // Annotating the tableau
   nonbasics_[0] = 0;
   for( size_t i=1; i<=n; ++i ) {
      basics_[i-1]  = -static_cast<ptrdiff_t>(i);
      nonbasics_[i] =  static_cast<ptrdiff_t>(i);
   }

   // Printing the annotated tableau
   //BLAZE_LOG_DEBUG_SECTION( log ) {
   //   printTableau( log );
   //}

   // Determination of the lexicographically smallest blocking variable for the initial pivot step
   size_t r = inf;

   for( size_t i=0; i<n; ++i ) {
      if( Q_(i,0) < -accuracy ) {  // < 0
         r = i;
         break;
      }
   }

   // We are finished if q >= 0 since z = 0 solves the LCP
   if( r == inf ) {
      lcp.x_ = real( 0 );
      return true;
   }

   for( size_t i=r+1; i<n; ++i ) {
      if( Q_(i,0) > -accuracy ) {  // >= 0
         BLAZE_INTERNAL_ASSERT( d[i] > -accuracy, "Negative value found" );  // >= 0
         continue;
      }

      BLAZE_INTERNAL_ASSERT( d[i] > real( 0 ), "Non-positive value found" );
      if( isLexicographicallyGreater( i, real(-1)/d[i], r, real(-1)/d[r]) )
         r = i;
   }

   size_t s = 0;
   size_t pivot_steps = 0;

   while( true )
   {
      // Perform pivot step
      pivot( r, s );

      // Printing the annotated tableau
      //BLAZE_LOG_DEBUG_SECTION( log ) {
      //   printTableau( log );
      //}

      // Finish if z0 blocked the driving variable
      if( nonbasics_[s] == 0 ) {
         lcp.x_ = real( 0 );
         for( size_t i=0; i<n; ++i ) {
            if( basics_[i] > 0 )
               lcp.x_[basics_[i]-1] = Q_(i,0);
         }

         return true;
      }

      // Find complement which will be driven next
      ptrdiff_t variable = -nonbasics_[s];
      s = inf;
      for( size_t i=0; i<=n; ++i ) {
         if( nonbasics_[i] == variable ) {
            s = i;
            break;
         }
      }

      BLAZE_INTERNAL_ASSERT( s != inf, "No complement found" );

      BLAZE_LOG_DEBUG_SECTION( log ) {
         if( nonbasics_[s] < 0 )
            log << "         w" << -nonbasics_[s] << " is new driving variable in column " << s;
         else
            log << "         z" <<  nonbasics_[s] << " is new driving variable in column " << s;
      }

      // Determination of the lexicographically smallest blocking variable
      r = inf;
      for( size_t i=0; i<n; ++i ) {
         if( M_(i,s) < -accuracy ) {  // < 0
            r = i;
            break;
         }
      }

      if( r == inf ) {
         // Driving variable is unblocked
         lcp.x_ = real( 0 );
         for( size_t i=0; i<n; ++i ) {
            if( basics_[i] > 0 )
               lcp.x_[basics_[i]-1] = Q_(i,0);
         }

         return false;
      }

      for( size_t i=r+1; i<n; ++i ) {
         if( M_(i,s) > -accuracy )  // >= 0
            continue;

         if( isLexicographicallyLess( i, real( -1 )/M_(i,s), r, real( -1 )/M_(r,s) ) )
            r = i;
      }

      ++pivot_steps;

      if( pivot_steps > 10*n ) {
         BLAZE_LOG_DEBUG_SECTION( log ) {
            log << BLAZE_YELLOW << "      WARNING: The LCP could not be solved with 10n pivot steps." << BLAZE_OLDCOLOR;
         }
         return false;
      }
   }

   return false;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief TODO
//
// \param v TODO
// \return TODO
//
// TODO < 0
*/
bool Lemke::isComponentwiseNonnegative( const VecN& v ) const
{
   for( size_t i=0; i<v.size(); ++i ) {
      if( v[i] < -accuracy )  // < 0
         return false;
   }
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
//
// \param v TODO
// \return TODO
//
// TODO <= 0
*/
bool Lemke::isComponentwisePositive( const VecN& v ) const
{
   for( size_t i=0; i<v.size(); ++i ) {
      if( v[i] < accuracy )  // <= 0
         return false;
   }
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
//
// \param i1 TODO
// \param f1 TODO
// \param i2 TODO
// \param f2 TODO
// \return TODO
//
// TODO
*/
bool Lemke::isLexicographicallyLess( size_t i1, real f1, size_t i2, real f2 ) const
{
   for( size_t j=0; j<Q_.columns(); ++j )
   {
      if( Q_(i1,j) * f1 < Q_(i2,j) * f2 - real(accuracy) )
         return true;
      else if( Q_(i1,j) * f1 > Q_(i2,j) * f2 + real(accuracy) )
         return false;
   }

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
//
// \param i1 TODO
// \param f1 TODO
// \param i2 TODO
// \param f2 TODO
// \return TODO
//
// TODO
*/
bool Lemke::isLexicographicallyGreater( size_t i1, real f1, size_t i2, real f2 ) const
{
   for( size_t j = 0; j < Q_.columns(); ++j )
   {
      if( Q_(i1,j) * f1 > Q_(i2,j) * f2 + real(accuracy) )
         return true;
      else if( Q_(i1,j) * f1 < Q_(i2,j) * f2 - real(accuracy) )
         return false;
   }

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
//
// \param block TODO
// \param drive TODO
// \return void
//
// TODO
*/
void Lemke::pivot( size_t block, size_t drive )
{
   using boost::format;

   size_t r( block );
   size_t s( drive );

   BLAZE_INTERNAL_ASSERT( r < basics_.size() && s <= basics_.size(), "Invalid preconditions for pivot step" );

   BLAZE_LOG_DEBUG_SECTION( log ) {
      log << "         pivot step <"
          << ( basics_[r] < 0 ? "w" : "z" )
          << ( basics_[r] < 0 ? -basics_[r] : basics_[r] ) << ", "
          << ( nonbasics_[s] < 0 ? "w" : "z" )
          << ( nonbasics_[s] < 0 ? -nonbasics_[s] : nonbasics_[s] ) << ">\n";
   }

   real invPivot = real( 1 ) / M_(r,s);

   for( size_t i=0; i<M_.rows(); ++i )
   {
      if( i == r ) continue;

      for( size_t j=0; j<Q_.columns(); ++j ) {
         Q_(i,j) -= Q_(r,j) * M_(i,s) * invPivot;
      }

      for( size_t j=0; j<M_.columns(); ++j )
      {
         if( j == s ) continue;
         M_(i,j) -= M_(i,s) * ( M_(r,j) * invPivot );
      }

      M_(i,s) *= invPivot;
   }

   for( size_t j=0; j<Q_.columns(); ++j ) {
      Q_(r,j) = -Q_(r,j) * invPivot;
   }

   for( size_t j=0; j<M_.columns(); ++j )
   {
      if( j == s )
         M_(r,j) = invPivot;
      else
         M_(r,j) = -M_(r,j) * invPivot;
   }

   // Swap the blocking and driving variables
   std::swap( basics_[r], nonbasics_[s] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
//
// \return TODO
//
// TODO
*/
void Lemke::printTableau( std::ostream& os ) const
{
   using boost::format;

   const size_t n( basics_.size() );

   // Printing the header
   os << "      ";
   for( size_t j=0; j<Q_.columns(); ++j ) {
      if( j == 0 )
         os << format( " %-7d " ) % 1;
      else
         os << format( " x%-6d " ) % j;
   }
   os << "  ";
   for( size_t j=0; j<M_.columns(); ++j ) {
      if( nonbasics_[j] < 0 )
         os << format( " w%-6d " ) % -nonbasics_[j];
      else
         os << format( " z%-6d " ) % nonbasics_[j];
   }
   os << "\n";

   // Printing the table border
   os << "     +-";
   for( size_t j=0; j<Q_.columns(); ++j ) {
      os << "---------";
   }
   os << "+-";
   for( size_t j=0; j<M_.columns(); ++j ) {
      os << "---------";
   }
   os << "+\n";

   for( size_t i=0; i<n; ++i ) {
      if( basics_[i] < 0 )
         os << format( " w%-2d | " ) % -basics_[i];
      else
         os << format( " z%-2d | " ) % basics_[i];

      for( size_t j=0; j<Q_.columns(); ++j ) {
         os << format( "%-8.2d " ) % Q_(i, j);
      }
      os << "| ";

      for( size_t j=0; j<M_.columns(); ++j ) {
         os << format( "%-8.2d " ) % M_(i, j);
      }
      os << "|\n";
   }

   os << "     +-";
   for( size_t j=0; j<Q_.columns(); ++j ) {
      os << "---------";
   }
   os << "+-";
   for( size_t j=0; j<M_.columns(); ++j ) {
      os << "---------";
   }
   os << "+\n\n";
}
//*************************************************************************************************

} // namespace blaze
