//=================================================================================================
/*!
//  \file blaze/math/MathTrait.h
//  \brief Header file for the mathematical/arithmetic trait
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

#ifndef _BLAZE_MATH_MATHTRAIT_H_
#define _BLAZE_MATH_MATHTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstddef>
#include <blaze/util/Complex.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>
#include <blaze/util/typetraits/RemoveCV.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  MATHEMATICAL TRAIT
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::MathTrait
// \brief Base template for the MathTrait class.
// \ingroup math
//
// \section mathtrait_general General
//
// The MathTrait class template offers the possibility to select the resulting data type of
// a generic mathematical operation between the two given types \a T1 and \a T2. In case of
// operations between built-in data types, the MathTrait class defines the more significant
// data type as the resulting data type. For this selection, larger and/or signed data types
// are given a higher significance. In case of operations involving user-defined data types,
// the MathTrait template specifies the resulting data type of this operation. \a const and
// \a volatile qualifiers and reference modifiers are generally ignored.\n
//
// MathTrait defines the following nested types:
//
//   - \a HighType : Represents the higher-order, more significant data type of the two
//                   given data types.
//   - \a LowType  : Represents the lower-order, less significant data type of the two
//                   given data types.
//   - \a AddType  : Represents the result of an addition operation.
//   - \a SubType  : Represents the result of a subtraction operation.
//   - \a MultType : Represents the result of a multiplication operation.
//   - \a CrossType: Represents the result of a cross product operation.
//   - \a DivType  : Represents the result of a division operation.
//
// If a certain mathematical operation is not possible and/or not defined between the two
// given data types, the according type is set to \a INVALID_TYPE.\n
//
// Specifying the resulting data type for a specific operation is done by specializing the
// MathTrait template for this particular type combination. In case a certain type combination
// is not defined in a MathTrait specialization, the base template is selected, which sets all
// nested types to \a INVALID_TYPE and therefore stops the compilation process. The following
// example shows the specialization for operations between the double and the integer type:

   \code
   template<>
   struct MathTrait< double, int >
   {
      typedef double        HighType;
      typedef int           LowType;
      typedef double        AddType;
      typedef double        SubType;
      typedef double        MultType;
      typedef INVALID_TYPE  CrossType;
      typedef double        DivType;
   };
   \endcode

// Per default, the MathTrait template provides specializations for the following built-in
// data types:
//
// <ul>
//    <li>integers</li>
//    <ul>
//       <li>unsigned char, signed char, char, wchar_t</li>
//       <li>unsigned short, short</li>
//       <li>unsigned int, int</li>
//       <li>unsigned long, long</li>
//       <li>std::size_t, std::ptrdiff_t (for certain 64-bit compilers)</li>
//    </ul>
//    <li>floating points</li>
//    <ul>
//       <li>float</li>
//       <li>double</li>
//       <li>long double</li>
//    </ul>
// </ul>
//
// Additionally, the Blaze library provides specializations for the following user-defined
// arithmetic types:
//
// <ul>
//    <li>std::complex</li>
//    <li>blaze::StaticVector</li>
//    <li>blaze::DynamicVector</li>
//    <li>blaze::CompressedVector</li>
//    <li>blaze::StaticMatrix</li>
//    <li>blaze::DynamicMatrix</li>
//    <li>blaze::CompressedMatrix</li>
//    <li>blaze::RotationMatrix</li>
//    <li>blaze::Quaternion</li>
// </ul>
//
//
// \n \section specializations Creating custom specializations
//
// It is possible to specialize the MathTrait template for additional user-defined data types.
// However, it is possible that a specific mathematical operation is invalid for the particular
// type combination. In this case, the \a INVALID_TYPE can be used to fill the missing type
// definition. The \a INVALID_TYPE represents the resulting data type of an invalid numerical
// operation. It is left undefined to stop the compilation process in case it is instantiated.
// The following example shows the specialization of the MathTrait template for 3D matrices and
// vectors. In this case, only the multiplication between the matrix and the vector is a valid
// numerical operation. Therefore for all other types the \a INVALID_TYPE is used.

   \code
   template< typename T1, typename T2 >
   struct MathTrait< StaticMatrix<T1,3UL,3UL,false>, StaticVector<T2,3UL,false> >
   {
      typedef typename typename MathTrait<T1,T2>::MultType  MT;

      typedef INVALID_TYPE                HighType;   // Invalid, no common high data type
      typedef INVALID_TYPE                LowType;    // Invalid, no common low data type
      typedef INVALID_TYPE                AddType;    // Invalid, cannot add a matrix and a vector
      typedef INVALID_TYPE                SubType;    // Invalid, cannot subtract a vector from a matrix
      typedef StaticVector<MT,3UL,false>  MultType;   // Multiplication between a matrix and a vector
      typedef INVALID_TYPE                CrossType;  // Invalid, cannot compute a matrix/vector cross product
      typedef INVALID_TYPE                DivType;    // Invalid, cannot divide a matrix by a vector
   };
   \endcode

// \n \section mathtrait_examples Examples
//
// The following example demonstrates the use of the MathTrait template, where depending on
// the two given data types the resulting data type is selected:

   \code
   template< typename T1, typename T2 >    // The two generic types
   typename MathTrait<T1,T2>::HighType     // The resulting generic return type
   add( T1 t1, T2 t2 )                     //
   {                                       // The function 'add' returns the sum
      return t1 + t2;                      // of the two given values
   }                                       //
   \endcode

// Additionally, the specializations of the MathTrait template enable arithmetic operations
// between any combination of the supported data types:

   \code
   // Vector of 3x3 single precision matrices
   typedef blaze::DynamicVector< blaze::StaticMatrix<float,3UL,3UL> >  VectorOfMatrices;

   // Vector of 3D double precision vectors
   typedef blaze::DynamicVector< blaze::StaticVector<double,3UL> >  VectorOfVectors;

   VectorOfMatrices vm;  // Setup of a vector of matrices
   VectorOfVectors  vv;  // Setup of a vector of vectors

   // Calculation of the component-wise vector product between the two vectors. The resulting
   // data type is again a dynamic vector of 3-dimensional vectors.
   VectorOfVectors res = vm * vv;
   \endcode
*/
template< typename T1, typename T2 >
struct MathTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure {
      typedef INVALID_TYPE  HighType;
      typedef INVALID_TYPE  LowType;
      typedef INVALID_TYPE  AddType;
      typedef INVALID_TYPE  SubType;
      typedef INVALID_TYPE  MultType;
      typedef INVALID_TYPE  CrossType;
      typedef INVALID_TYPE  DivType;
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { qualified = IsConst<T1>::value || IsVolatile<T1>::value || IsReference<T1>::value ||
                      IsConst<T2>::value || IsVolatile<T2>::value || IsReference<T2>::value };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename RemoveReference< typename RemoveCV<T1>::Type >::Type  Type1;
   typedef typename RemoveReference< typename RemoveCV<T2>::Type >::Type  Type2;
   typedef MathTrait<Type1,Type2>  Helper;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, Helper, Failure >::Type::HighType   HighType;
   typedef typename SelectType< qualified, Helper, Failure >::Type::LowType    LowType;
   typedef typename SelectType< qualified, Helper, Failure >::Type::AddType    AddType;
   typedef typename SelectType< qualified, Helper, Failure >::Type::SubType    SubType;
   typedef typename SelectType< qualified, Helper, Failure >::Type::MultType   MultType;
   typedef typename SelectType< qualified, Helper, Failure >::Type::CrossType  CrossType;
   typedef typename SelectType< qualified, Helper, Failure >::Type::DivType    DivType;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  MATHTRAIT SPECIALIZATION MACRO
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the creation of MathTrait specializations for the built-in data types.
// \ingroup math
//
// This macro is used for the setup of the MathTrait specializations for the built-in data types.
*/
#define BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION(T1,T2,HIGH,LOW) \
   template<> \
   struct MathTrait< T1, T2 > \
   { \
      typedef HIGH          HighType; \
      typedef LOW           LowType;  \
      typedef HIGH          AddType;  \
      typedef HIGH          SubType;  \
      typedef HIGH          MultType; \
      typedef INVALID_TYPE  CrossType; \
      typedef HIGH          DivType;  \
   }
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the creation of MathTrait specializations for the complex data type.
// \ingroup math
//
// This macro is used for the setup of the MathTrait specializations for the complex data type.
*/
#define BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( T1 ) \
   template< typename T2 > \
   struct MathTrait< T1, complex<T2> > \
   { \
      typedef complex<T2>   HighType; \
      typedef T1            LowType;  \
      typedef complex<T2>   AddType;  \
      typedef complex<T2>   SubType;  \
      typedef complex<T2>   MultType; \
      typedef INVALID_TYPE  CrossType; \
      typedef complex<T2>   DivType;  \
   }; \
   template< typename T2 > \
   struct MathTrait< complex<T2>, T1 > \
   { \
      typedef complex<T2>   HighType; \
      typedef T1            LowType;  \
      typedef complex<T2>   AddType;  \
      typedef complex<T2>   SubType;  \
      typedef complex<T2>   MultType; \
      typedef INVALID_TYPE  CrossType; \
      typedef complex<T2>   DivType;  \
   }
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UNSIGNED CHAR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , unsigned char , unsigned char , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , char          , char          , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , signed char   , signed char   , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , wchar_t       , wchar_t       , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , unsigned short, unsigned short, unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , short         , short         , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , unsigned int  , unsigned int  , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , int           , int           , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , unsigned long , unsigned long , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , long          , long          , unsigned char  );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , std::size_t   , std::size_t   , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , std::ptrdiff_t, std::ptrdiff_t, unsigned char  );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , float         , float         , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , double        , double        , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned char , long double   , long double   , unsigned char  );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CHAR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , unsigned char , char          , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , char          , char          , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , signed char   , signed char   , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , wchar_t       , wchar_t       , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , unsigned short, unsigned short, char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , short         , short         , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , unsigned int  , unsigned int  , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , int           , int           , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , unsigned long , unsigned long , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , long          , long          , char           );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , std::size_t   , std::size_t   , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , std::ptrdiff_t, std::ptrdiff_t, char           );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , float         , float         , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , double        , double        , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( char          , long double   , long double   , char           );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIGNED CHAR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , unsigned char , signed char   , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , char          , signed char   , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , signed char   , signed char   , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , wchar_t       , wchar_t       , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , unsigned short, unsigned short, signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , short         , short         , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , unsigned int  , unsigned int  , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , int           , int           , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , unsigned long , unsigned long , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , long          , long          , signed char    );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , std::size_t   , std::size_t   , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , std::ptrdiff_t, std::ptrdiff_t, signed char    );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , float         , float         , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , double        , double        , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( signed char   , long double   , long double   , signed char    );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  WCHAR_T SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , unsigned char , wchar_t       , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , char          , wchar_t       , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , signed char   , wchar_t       , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , wchar_t       , wchar_t       , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , unsigned short, unsigned short, wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , short         , short         , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , unsigned int  , unsigned int  , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , int           , int           , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , unsigned long , unsigned long , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , long          , long          , wchar_t        );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , std::size_t   , std::size_t   , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , std::ptrdiff_t, std::ptrdiff_t, wchar_t        );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , float         , float         , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , double        , double        , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( wchar_t       , long double   , long double   , wchar_t        );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UNSIGNED SHORT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, unsigned char , unsigned short, unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, char          , unsigned short, char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, signed char   , unsigned short, signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, wchar_t       , unsigned short, wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, unsigned short, unsigned short, unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, short         , short         , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, unsigned int  , unsigned int  , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, int           , int           , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, unsigned long , unsigned long , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, long          , long          , unsigned short );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, std::size_t   , std::size_t   , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, std::ptrdiff_t, std::ptrdiff_t, unsigned short );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, float         , float         , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, double        , double        , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned short, long double   , long double   , unsigned short );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SHORT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , unsigned char , short         , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , char          , short         , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , signed char   , short         , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , wchar_t       , short         , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , unsigned short, short         , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , short         , short         , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , unsigned int  , unsigned int  , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , int           , int           , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , unsigned long , unsigned long , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , long          , long          , short          );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , std::size_t   , std::size_t   , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , std::ptrdiff_t, std::ptrdiff_t, short          );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , float         , float         , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , double        , double        , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( short         , long double   , long double   , short          );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UNSIGNED INT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , unsigned char , unsigned int  , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , char          , unsigned int  , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , signed char   , unsigned int  , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , wchar_t       , unsigned int  , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , unsigned short, unsigned int  , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , short         , unsigned int  , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , unsigned int  , unsigned int  , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , int           , int           , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , unsigned long , unsigned long , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , long          , long          , unsigned int   );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , std::size_t   , std::size_t   , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , std::ptrdiff_t, std::ptrdiff_t, unsigned int   );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , float         , float         , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , double        , double        , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned int  , long double   , long double   , unsigned int   );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , unsigned char , int           , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , char          , int           , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , signed char   , int           , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , wchar_t       , int           , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , unsigned short, int           , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , short         , int           , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , unsigned int  , int           , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , int           , int           , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , unsigned long , unsigned long , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , long          , long          , int            );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , std::size_t   , std::size_t   , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , std::ptrdiff_t, std::ptrdiff_t, int            );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , float         , float         , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , double        , double        , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( int           , long double   , long double   , int            );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UNSIGNED LONG SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , unsigned char , unsigned long , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , char          , unsigned long , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , signed char   , unsigned long , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , wchar_t       , unsigned long , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , unsigned short, unsigned long , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , short         , unsigned long , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , unsigned int  , unsigned long , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , int           , unsigned long , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , unsigned long , unsigned long , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , long          , long          , unsigned long  );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , std::size_t   , std::size_t   , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , std::ptrdiff_t, std::ptrdiff_t, unsigned long  );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , float         , float         , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , double        , double        , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( unsigned long , long double   , long double   , unsigned long  );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LONG SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , unsigned char , long          , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , char          , long          , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , signed char   , long          , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , wchar_t       , long          , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , unsigned short, long          , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , short         , long          , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , unsigned int  , long          , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , int           , long          , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , unsigned long , long          , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , long          , long          , long           );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , std::size_t   , std::size_t   , long           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , std::ptrdiff_t, std::ptrdiff_t, long           );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , float         , float         , long           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , double        , double        , long           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long          , long double   , long double   , long           );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE_T SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
#if defined(_WIN64)
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , unsigned char , std::size_t   , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , char          , std::size_t   , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , signed char   , std::size_t   , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , wchar_t       , std::size_t   , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , unsigned short, std::size_t   , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , short         , std::size_t   , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , unsigned int  , std::size_t   , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , int           , std::size_t   , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , unsigned long , std::size_t   , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , long          , std::size_t   , long           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , std::size_t   , std::size_t   , std::size_t    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , std::ptrdiff_t, std::ptrdiff_t, std::size_t    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , float         , float         , std::size_t    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , double        , double        , std::size_t    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::size_t   , long double   , long double   , std::size_t    );
/*! \endcond */
#endif
//*************************************************************************************************




//=================================================================================================
//
//  PTRDIFF_T SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
#if defined(_WIN64)
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, unsigned char , std::ptrdiff_t, unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, char          , std::ptrdiff_t, char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, signed char   , std::ptrdiff_t, signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, wchar_t       , std::ptrdiff_t, wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, unsigned short, std::ptrdiff_t, unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, short         , std::ptrdiff_t, short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, unsigned int  , std::ptrdiff_t, unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, int           , std::ptrdiff_t, int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, unsigned long , std::ptrdiff_t, unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, long          , std::ptrdiff_t, long           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, std::size_t   , std::ptrdiff_t, std::size_t    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, float         , float         , std::ptrdiff_t );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, double        , double        , std::ptrdiff_t );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t, long double   , long double   , std::ptrdiff_t );
/*! \endcond */
#endif
//*************************************************************************************************




//=================================================================================================
//
//  FLOAT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , unsigned char , float         , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , char          , float         , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , signed char   , float         , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , wchar_t       , float         , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , unsigned short, float         , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , short         , float         , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , unsigned int  , float         , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , int           , float         , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , unsigned long , float         , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , long          , float         , long           );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , std::size_t   , float         , std::size_t    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , std::ptrdiff_t, float         , std::ptrdiff_t );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , float         , float         , float          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , double        , double        , float          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( float         , long double   , long double   , float          );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DOUBLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , unsigned char , double        , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , char          , double        , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , signed char   , double        , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , wchar_t       , double        , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , unsigned short, double        , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , short         , double        , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , unsigned int  , double        , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , int           , double        , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , unsigned long , double        , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , long          , double        , long           );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , std::size_t   , double        , std::size_t    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , std::ptrdiff_t, double        , std::ptrdiff_t );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , float         , double        , float          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , double        , double        , double         );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( double        , long double   , long double   , double         );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LONG DOUBLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//                                             Type 1          Type 2          High type       Low type
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , unsigned char , long double   , unsigned char  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , char          , long double   , char           );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , signed char   , long double   , signed char    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , wchar_t       , long double   , wchar_t        );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , unsigned short, long double   , unsigned short );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , short         , long double   , short          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , unsigned int  , long double   , unsigned int   );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , int           , long double   , int            );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , unsigned long , long double   , unsigned long  );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , long          , long double   , long           );
#if defined(_WIN64)
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , std::size_t   , long double   , std::size_t    );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , std::ptrdiff_t, long double   , std::ptrdiff_t );
#endif
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , float         , long double   , float          );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , double        , long double   , double         );
BLAZE_CREATE_BUILTIN_MATHTRAIT_SPECIALIZATION( long double   , long double   , long double   , long double    );
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COMPLEX SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( unsigned char  );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( char           );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( signed char    );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( wchar_t        );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( unsigned short );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( short          );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( unsigned int   );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( int            );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( unsigned long  );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( long           );
#if defined(_WIN64)
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( std::size_t    );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( std::ptrdiff_t );
#endif
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( float          );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( double         );
BLAZE_CREATE_COMPLEX_MATHTRAIT_SPECIALIZATION( long double    );
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct MathTrait< complex<T1>, complex<T2> >
{
   typedef complex<typename MathTrait<T1,T2>::HighType>  HighType;
   typedef complex<typename MathTrait<T1,T2>::LowType>   LowType;
   typedef complex<typename MathTrait<T1,T2>::AddType>   AddType;
   typedef complex<typename MathTrait<T1,T2>::SubType>   SubType;
   typedef complex<typename MathTrait<T1,T2>::MultType>  MultType;
   typedef complex<typename MathTrait<T1,T2>::DivType>   DivType;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
