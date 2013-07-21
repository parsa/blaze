//=================================================================================================
/*!
//  \file blaze/math/serialization/VectorSerializer.h
//  \brief Serialization of dense and sparse vectors
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

#ifndef _BLAZE_MATH_SERIALIZATION_VECTORSERIALIZER_H_
#define _BLAZE_MATH_SERIALIZATION_VECTORSERIALIZER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <blaze/math/constraints/Vector.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/Vector.h>
#include <blaze/math/serialization/TypeValueMapping.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Serializer for dense and sparse vectors.
// \ingroup math_serialization
//
// The VectorSerializer implements the necessary logic to serialize dense and sparse vectors, i.e.
// to convert them into a portable, binary representation. The following example demonstrates the
// (de-)serialization process of vectors:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   // Serialization of both vectors
   {
      blaze::StaticVector<double,5UL,rowVector> d;
      blaze::CompressedVector<int,columnVector> s;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "vectors.blaze"
      blaze::Archive<std::ofstream> archive( "vectors.blaze" );

      // Serialization of both vectors into the same archive. Note that d lies before s!
      archive << d << s;
   }

   // Reconstitution of both vectors
   {
      blaze::DynamicVector<double,rowVector> d1;
      blaze::DynamicVector<int,rowVector> d2;

      // Creating an archive that reads from the file "vectors.blaze"
      blaze::Archive<std::ofstream> archive( "vectors.blaze" );

      // Reconstituting the former d vector into d1. Note that it is possible to reconstitute
      // the vector into a differrent kind of vector (StaticVector -> DynamicVector), but that
      // the type of elements has to be the same.
      archive >> d1;

      // Reconstituting the former s vector into d2. Note that is is even possible to reconstitute
      // a sparse vector as a dense vector (also the reverse is possible) and that a column vector
      // can be reconstituted as row vector (and vice versa). Note however that also in this case
      // the type of elements is the same!
      archive >> d2
   }
   \endcode

// Note that it is even possible to (de-)serialize vectors with vector or matrix elements:

   \code
   // Serialization
   {
      blaze::CompressedVector< blaze::DynamicVector< blaze::complex<double> > > vec;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "vector.blaze"
      blaze::Archive<std::ofstream> archive( "vector.blaze" );

      // Serialization of the vector into the archive
      archive << vec;
   }

   // Deserialization
   {
      blaze::CompressedVector< blaze::DynamicVector< blaze::complex<double> > > vec;

      // Creating an archive that reads from the file "vector.blaze"
      blaze::Archive<std::ofstream> archive( "vector.blaze" );

      // Reconstitution of the vector from the archive
      archive >> vec;
   }
   \endcode

// As the examples demonstrates, the vector serialization offers an enormous flexibility. However,
// several actions result in errors:
//
//  - vectors cannot be reconstituted as matrices (and vice versa)
//  - the element type of the serialized and reconstituted vector must match, which means
//    that on the source and destination platform the general type (signed/unsigned integral
//    or floating point) and the size of the type must be exactly the same
//  - when reconstituting a StaticVector, its size must match the size of the serialized vector
//
// In case an error is encountered during (de-)serialization, a \a std::runtime_exception is
// thrown.
*/
class VectorSerializer
{
 private:
   //**Private class VectorValueMappingHelper******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Auxiliary helper class for the VectorValueMapping class template.
   */
   template< bool IsDenseVector >
   struct VectorValueMappingHelper;
   /*! \endcond */
   //**********************************************************************************************

   //**Private class VectorValueMapping************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Serialization of the type of a vector.
   //
   // This class template converts the given vector type into an integral representation suited
   // for serialization. Depending on the given vector type, the \a value member enumeration is
   // set to the according integral representation.
   */
   template< typename T >
   struct VectorValueMapping
   {
      enum { value = VectorValueMappingHelper< IsDenseVector<T>::value >::value };
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE( T );
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline VectorSerializer();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   // No explicitly declared copy assignment operator.
   //**********************************************************************************************

   //**Serialization functions*********************************************************************
   /*!\name Serialization functions */
   //@{
   template< typename Archive, typename VT, bool TF >
   void serialize( Archive& archive, const Vector<VT,TF>& vec );
   //@}
   //**********************************************************************************************

   //**Deserialization functions*********************************************************************
   /*!\name Deserialization functions */
   //@{
   template< typename Archive, typename VT, bool TF >
   void deserialize( Archive& archive, Vector<VT,TF>& vec );
   //@}
   //**********************************************************************************************

 private:
   //**Serialization functions*********************************************************************
   /*!\name Serialization functions */
   //@{
   template< typename Archive, typename VT >
   void serializeHeader( Archive& archive, const VT& vec );

   template< typename Archive, typename VT, bool TF >
   void serializeVector( Archive& archive, const DenseVector<VT,TF>& vec );

   template< typename Archive, typename VT, bool TF >
   void serializeVector( Archive& archive, const SparseVector<VT,TF>& vec );
   //@}
   //**********************************************************************************************

   //**Deserialization functions*******************************************************************
   /*!\name Deserialization functions */
   //@{
   template< typename Archive, typename VT >
   void deserializeHeader( Archive& archive, const VT& vec );

   template< typename VT >
   typename DisableIf< IsResizable<VT> >::Type prepareVector( VT& vec );

   template< typename VT >
   typename EnableIf< IsResizable<VT> >::Type prepareVector( VT& vec );

   template< typename Archive, typename VT >
   void deserializeVector( Archive& archive, VT& vec );

   template< typename Archive, typename VT, bool TF >
   typename DisableIfTrue< IsNumeric< typename VT::ElementType >::value && VT::vectorizable >::Type
      deserializeDenseVector( Archive& archive, DenseVector<VT,TF>& vec );

   template< typename Archive, typename VT, bool TF >
   typename EnableIfTrue< IsNumeric< typename VT::ElementType >::value && VT::vectorizable >::Type
      deserializeDenseVector( Archive& archive, DenseVector<VT,TF>& vec );

   template< typename Archive, typename VT, bool TF >
   void deserializeDenseVector( Archive& archive, SparseVector<VT,TF>& vec );

   template< typename Archive, typename VT, bool TF >
   void deserializeSparseVector( Archive& archive, DenseVector<VT,TF>& vec );

   template< typename Archive, typename VT, bool TF >
   void deserializeSparseVector( Archive& archive, SparseVector<VT,TF>& vec );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   uint8_t  version_;      //!< The version of the archive.
   uint8_t  type_;         //!< The type of the vector.
   uint8_t  elementType_;  //!< The type of an element.
   uint8_t  elementSize_;  //!< The size in bytes of a single element of the vector.
   uint64_t size_;         //!< The size of the vector.
   uint64_t number_;       //!< The total number of elements contained in the vector.
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
/*!\brief The default constructor of the VectorSerializer class.
*/
VectorSerializer::VectorSerializer()
   : version_    ( 0U  )  // The version of the archive
   , type_       ( 0U  )  // The type of the vector
   , elementType_( 0U  )  // The type of an element
   , elementSize_( 0U  )  // The size in bytes of a single element of the vector
   , size_       ( 0UL )  // The size of the vector
   , number_     ( 0UL )  // The total number of elements contained in the vector
{}
//*************************************************************************************************




//=================================================================================================
//
//  SERIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Serializes the given vector and writes it to the archive.
//
// \param archive The archive to be written.
// \param vec The vector to be serialized.
// \return void
// \exception std::runtime_error Error during serialization.
//
// This function serializes the given vector and writes it to the given archive. In case any
// error is detected during the serialization, a \a std::runtime_error is thrown.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void VectorSerializer::serialize( Archive& archive, const Vector<VT,TF>& vec )
{
   if( !archive ) {
      throw std::runtime_error( "Faulty archive detected" );
   }

   serializeHeader( archive, ~vec );
   serializeVector( archive, ~vec );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Serializes all meta information about the given vector.
//
// \param archive The archive to be written.
// \param vec The vector to be serialized.
// \return void
// \exception std::runtime_error File header could not be serialized.
*/
template< typename Archive  // Type of the archive
        , typename VT >     // Type of the vector
void VectorSerializer::serializeHeader( Archive& archive, const VT& vec )
{
   typedef typename VT::ElementType  ET;

   archive << uint8_t ( 1U );
   archive << uint8_t ( VectorValueMapping<VT>::value );
   archive << uint8_t ( TypeValueMapping<ET>::value );
   archive << uint8_t ( sizeof( ET ) );
   archive << uint64_t( vec.size() );
   archive << uint64_t( IsDenseVector<VT>::value ? vec.size() : vec.nonZeros() );

   if( !archive ) {
      throw std::runtime_error( "File header could not be serialized" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Serializes the elements of a dense vector.
//
// \param archive The archive to be written.
// \param vec The vector to be serialized.
// \return void
// \exception std::runtime_error Dense vector could not be serialized.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void VectorSerializer::serializeVector( Archive& archive, const DenseVector<VT,TF>& vec )
{
   size_t i( 0UL );
   while( ( i < (~vec).size() ) && ( archive << (~vec)[i] ) ) {
      ++i;
   }

   if( !archive ) {
      throw std::runtime_error( "Dense vector could not be serialized" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Serializes the elements of a sparse vector.
//
// \param archive The archive to be written.
// \param vec The vector to be serialized.
// \return void
// \exception std::runtime_error Sparse vector could not be serialized.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void VectorSerializer::serializeVector( Archive& archive, const SparseVector<VT,TF>& vec )
{
   typedef typename VT::ConstIterator  ConstIterator;

   ConstIterator element( (~vec).begin() );
   while( ( element != (~vec).end() ) &&
          ( archive << element->index() << element->value() ) ) {
      ++element;
   }

   if( !archive ) {
      throw std::runtime_error( "Sparse vector could not be serialized" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DESERIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Deserializes a vector from the given archive.
//
// \param archive The archive to be read from.
// \param vec The vector to be deserialized.
// \return void
// \exception std::runtime_error Error during deserialization.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void VectorSerializer::deserialize( Archive& archive, Vector<VT,TF>& vec )
{
   if( !archive ) {
      throw std::invalid_argument( "Faulty archive detected" );
   }

   deserializeHeader( archive, ~vec );
   prepareVector( ~vec );
   deserializeVector( archive, ~vec );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deserializes all meta information about the given vector.
//
// \param archive The archive to be read from.
// \param vec The vector to be deserialized.
// \return void
// \exception std::runtime_error Error during deserialization.
//
// This function deserializes all meta information about the given vector contained in the
// header of the given archive. In case any error is detected during the deserialization
// process (for instance an invalid type of vector, element type, element size, or vector
// size) a \a std::runtime_error is thrown.
*/
template< typename Archive  // Type of the archive
        , typename VT >     // Type of the vector
void VectorSerializer::deserializeHeader( Archive& archive, const VT& vec )
{
   typedef typename VT::ElementType  ET;

   if( !( archive >> version_ >> type_ >> elementType_ >> elementSize_ >> size_ >> number_ ) ) {
      throw std::runtime_error( "Corrupt archive detected" );
   }
   else if( version_ != 1UL ) {
      throw std::runtime_error( "Invalid version detected" );
   }
   else if( type_ != 0U && type_ != 1U ) {
      throw std::runtime_error( "Invalid vector type detected" );
   }
   else if( elementType_ != TypeValueMapping<ET>::value ) {
      throw std::runtime_error( "Invalid element type detected" );
   }
   else if( elementSize_ != sizeof( ET ) ) {
      throw std::runtime_error( "Invalid element size detected" );
   }
   else if( !IsResizable<VT>::value && size_ != vec.size() ) {
      throw std::runtime_error( "Invalid vector size detected" );
   }
   else if( number_ > size_ ) {
      throw std::runtime_error( "Invalid number of elements detected" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Prepares the given non-resizable vector for the deserialization process.
//
// \param vec The vector to be prepared.
// \return void
*/
template< typename VT >  // Type of the vector
typename DisableIf< IsResizable<VT> >::Type VectorSerializer::prepareVector( VT& vec )
{
   reset( vec );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Prepares the given resizable vector for the deserialization process.
//
// \param vec The vector to be prepared.
// \return void
*/
template< typename VT >  // Type of the vector
typename EnableIf< IsResizable<VT> >::Type VectorSerializer::prepareVector( VT& vec )
{
   vec.resize ( size_, false );
   vec.reserve( number_ );
   reset( vec );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deserializes a vector from the archive.
//
// \param archive The archive to be read from.
// \param vec The vector to be reconstituted.
// \return void
// \exception std::runtime_error Error during deserialization.
//
// This function deserializes the contents of the vector from the archive and reconstitutes the
// given vector.
*/
template< typename Archive  // Type of the archive
        , typename VT >     // Type of the vector
void VectorSerializer::deserializeVector( Archive& archive, VT& vec )
{
   if( type_ == 0U ) {
      deserializeDenseVector( archive, vec );
   }
   else if( type_ == 1U ) {
      deserializeSparseVector( archive, vec );
   }
   else {
      BLAZE_INTERNAL_ASSERT( false, "Undefined type flag" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deserializes a dense vector from the archive.
//
// \param archive The archive to be read from.
// \param vec The dense vector to be reconstituted.
// \return void
// \exception std::runtime_error Dense vector could not be deserialized.
//
// This function deserializes a dense vector from the archive and reconstitutes the given
// dense vector. In case any error is detected during the deserialization process, a
// \a std::runtime_error is thrown.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
typename DisableIfTrue< IsNumeric< typename VT::ElementType >::value && VT::vectorizable >::Type
   VectorSerializer::deserializeDenseVector( Archive& archive, DenseVector<VT,TF>& vec )
{
   typedef typename VT::ElementType  ET;

   size_t i( 0UL );
   ET value = ET();

   while( ( i != size_ ) && ( archive >> value ) ) {
      (~vec)[i] = value;
      ++i;
   }

   if( !archive ) {
      throw std::runtime_error( "Dense vector could not be deserialized" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deserializes a dense vector from the archive.
//
// \param archive The archive to be read from.
// \param vec The dense vector to be reconstituted.
// \return void
// \exception std::runtime_error Dense vector could not be deserialized.
//
// This function deserializes a dense vector from the archive and reconstitutes the given
// dense vector. In case any error is detected during the deserialization process, a
// \a std::runtime_error is thrown.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
typename EnableIfTrue< IsNumeric< typename VT::ElementType >::value && VT::vectorizable >::Type
   VectorSerializer::deserializeDenseVector( Archive& archive, DenseVector<VT,TF>& vec )
{
   if( size_ == 0UL ) return;
   archive.read( &(~vec)[0], size_ );

   if( !archive ) {
      throw std::runtime_error( "Dense vector could not be deserialized" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deserializes a dense vector from the archive.
//
// \param archive The archive to be read from.
// \param vec The sparse vector to be reconstituted.
// \return void
// \exception std::runtime_error Sparse vector could not be deserialized.
//
// This function deserializes a dense vector from the archive and reconstitutes the given
// sparse vector. In case any error is detected during the deserialization process, a
// \a std::runtime_error is thrown.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void VectorSerializer::deserializeDenseVector( Archive& archive, SparseVector<VT,TF>& vec )
{
   typedef typename VT::ElementType  ET;

   size_t i( 0UL );
   ET value = ET();

   while( ( i != size_ ) && ( archive >> value ) ) {
      (~vec)[i] = value;
      ++i;
   }

   if( !archive ) {
      throw std::runtime_error( "Sparse vector could not be deserialized" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deserializes a sparse vector from the archive.
//
// \param archive The archive to be read from.
// \param vec The dense vector to be reconstituted.
// \return void
// \exception std::runtime_error Dense vector could not be deserialized.
//
// This function deserializes a sparse vector from the archive and reconstitutes the given
// dense vector. In case any error is detected during the deserialization process, a
// \a std::runtime_error is thrown.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void VectorSerializer::deserializeSparseVector( Archive& archive, DenseVector<VT,TF>& vec )
{
   typedef typename VT::ElementType  ET;

   size_t i( 0UL );
   size_t index( 0UL );
   ET     value = ET();

   while( ( i != number_ ) && ( archive >> index >> value ) ) {
      (~vec)[index] = value;
      ++i;
   }

   if( !archive ) {
      throw std::runtime_error( "Dense vector could not be deserialized" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deserializes a sparse vector from the archive.
//
// \param archive The archive to be read from.
// \param vec The sparse vector to be reconstituted.
// \return void
// \exception std::runtime_error Sparse vector could not be deserialized.
//
// This function deserializes a sparse vector from the archive and reconstitutes the given
// sparse vector. In case any error is detected during the deserialization process, a
// \a std::runtime_error is thrown.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void VectorSerializer::deserializeSparseVector( Archive& archive, SparseVector<VT,TF>& vec )
{
   typedef typename VT::ElementType  ET;

   size_t i( 0UL );
   size_t index( 0UL );
   ET     value = ET();

   while( ( i != number_ ) && ( archive >> index >> value ) ) {
      (~vec).append( index, value, false );
      ++i;
   }

   if( !archive ) {
      throw std::runtime_error( "Sparse vector could not be deserialized" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  VECTORVALUEMAPPINGHELPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the VectorValueMappingHelper class template for dense vectors.
*/
template<>
struct VectorSerializer::VectorValueMappingHelper<true>
{
   enum { value = 0 };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the VectorValueMappingHelper class template for sparse vectors.
*/
template<>
struct VectorSerializer::VectorValueMappingHelper<false>
{
   enum { value = 1 };
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Serializes the given vector and writes it to the archive.
//
// \param archive The archive to be written.
// \param vec The vector to be serialized.
// \return void
// \exception std::runtime_error Error during serialization.
//
// The serialize() function converts the given vector into a portable, binary representation.
// The following example demonstrates the (de-)serialization process of vectors:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   // Serialization of both vectors
   {
      blaze::StaticVector<double,5UL,rowVector> d;
      blaze::CompressedVector<int,columnVector> s;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "vectors.blaze"
      blaze::Archive<std::ofstream> archive( "vectors.blaze" );

      // Serialization of both vectors into the same archive. Note that d lies before s!
      archive << d << s;
   }

   // Reconstitution of both vectors
   {
      blaze::DynamicVector<double,rowVector> d1;
      blaze::DynamicVector<int,rowVector> d2;

      // ... Resizing and initialization

      // Creating an archive that reads from the file "vectors.blaze"
      blaze::Archive<std::ofstream> archive( "vectors.blaze" );

      // Reconstituting the former d vector into d1. Note that it is possible to reconstitute
      // the vector into a differrent kind of vector (StaticVector -> DynamicVector), but that
      // the type of elements has to be the same.
      archive >> d1;

      // Reconstituting the former s vector into d2. Note that is is even possible to reconstitute
      // a sparse vector as a dense vector (also the reverse is possible) and that a column vector
      // can be reconstituted as row vector (and vice versa). Note however that also in this case
      // the type of elements is the same!
      archive >> d2
   }
   \endcode

// As the example demonstrates, the vector serialization offers an enormous flexibility. However,
// several actions result in errors:
//
//  - vectors cannot be reconstituted as matrices (and vice versa)
//  - the element type of the serialized and reconstituted vector must match, which means
//    that on the source and destination platform the general type (signed/unsigned integral
//    or floating point) and the size of the type must be exactly the same
//  - when reconstituting a StaticVector, its size must match the size of the serialized vector
//
// In case an error is encountered during (de-)serialization, a \a std::runtime_exception is
// thrown.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void serialize( Archive& archive, const Vector<VT,TF>& vec )
{
   VectorSerializer().serialize( archive, ~vec );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deserializes a vector from the given archive.
//
// \param archive The archive to be read from.
// \param vec The vector to be deserialized.
// \return void
// \exception std::runtime_error Vector could not be deserialized.
//
// The deserialize() function converts the portable, binary representation contained in
// the given archive into the given vector type. For a detailed example that demonstrates
// the (de-)serialization process of vectors, see the serialize() function.
*/
template< typename Archive  // Type of the archive
        , typename VT       // Type of the vector
        , bool TF >         // Transpose flag
void deserialize( Archive& archive, Vector<VT,TF>& vec )
{
   VectorSerializer().deserialize( archive, ~vec );
}
//*************************************************************************************************

} // namespace blaze

#endif
