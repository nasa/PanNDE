/*! \headerfile HostArray.hpp "modules/HostData/include/HostArray.hpp"
* "HostArray.hpp" contains class implementations for host-based array management,
* data bundling, and array factory functionality. These classes provide standard 
* interfaces for array operations while avoiding dependencies on the STL in public
* interfaces to maintain compatibility with device-based extensions (like GPU code).
*
*/

/******************************************************
# Notices:
# Copyright 2021 United States Government as represented by the 
# Administrator of the National Aeronautics and Space Administration. 
# No copyright is claimed in the United States under Title 17, U.S. Code. 
# All Other Rights Reserved.
# 
# Googletest is a product of Google Inc. and is subject to the following:
# 
# Copyright 2008, Google Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice, 
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
#  * Neither the name of Google Inc. nor the names of its contributors may be used 
#    to endorse or promote products derived from this software without specific 
#    prior written permission.
#
# Disclaimers
# No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, 
# EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT 
# THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT 
# SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO 
# THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY 
# GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE 
# PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT 
# AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE 
# ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS." 
#
# Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES 
# GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S 
# USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES 
# ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S 
# USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, 
# ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  
# RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS 
# AGREEMENT.
******************************************************/

#pragma once

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include <vector>
#include <map>

#include "Array.hpp"

namespace HostData {
  /*! \class HostArray<T> HostArray.hpp "modules/HostData/include/HostArray.hpp"
  *
  * Implements a fixed-size array class providing the PanNDE::Array interface.
  * This class enables cross-platform compatibility by avoiding STL dependencies
  * in public interfaces. Internally it uses std::vector for storage while providing
  * a consistent interface that can be implemented for various hardware targets.
  *
  */
  template<typename T>
  class HostArray : public PanNDE::Array<T> {
    public:
      /*!
      * Creates a shared pointer to a new HostArray instance.
      * \return std::shared_ptr<HostData::HostArray<T>> Shared pointer to the newly created array
      */
      static
      std::shared_ptr<HostData::HostArray<T>> makeShared(){
        auto array=std::make_shared<HostData::HostArray<T>>(HostData::HostArray<T>());
        return std::move(array);
      };
      
      /*!
      * Accesses element at specified index with bounds checking.
      * \param idx int64_t Index of the element to access
      * \return T& Reference to the element at the specified position
      * \throw std::out_of_range If idx is out of bounds of the array
      */
      T& at(int64_t idx)override{return vec.at(idx);};
      
      /*!
      * Returns the number of elements in the array.
      * \return int64_t The number of elements in the array
      */
      int64_t size()override{return vec.size();};
      
      /*!
      * Resizes the array to contain the specified number of elements.
      * Note: Array contents are not guaranteed to be preserved after resizing.
      * \param newsize int64_t New size of the array
      */
      void resize(int64_t newsize)override{return vec.resize(newsize);};
      
      /*!
      * Provides direct access to the underlying data buffer.
      * Only use when low-level access is required by performance-critical code.
      * \return T* Pointer to the underlying data array
      */
      T* data()override{return vec.data();};

      /*!
      * Copies data from another array into this array.
      * \param other PanNDE::Array<T>* Pointer to source array to copy from
      */
      void copy(PanNDE::Array<T>* other)override{
        this->resize(other->size());
        for(auto k=0;k<(this->size());k++){this->at(k)=other->at(k);};
      };
      
      /*!
      * Copies data from another array (shared pointer) into this array.
      * \param other std::shared_ptr<PanNDE::Array<T>> Shared pointer to source array
      */
      void copy(std::shared_ptr<PanNDE::Array<T>> other)override{
        copy(other.get());
      };

      /*!
      * Copies data from a C-style array structure into this array.
      * This method is primarily intended for inter-process or device communications.
      * \param other PanNDE::CArray<T>& Reference to C-array structure to copy from
      * \warning CArray has no bounds checking, use with caution
      */
      void copy(PanNDE::CArray<T>& other)override{
        this->resize(other.length);
        for(auto k=0;k<(this->size());k++){this->at(k)=other.data[k];};
      };
      
      /*!
      * Returns a C-style array representation of the array data.
      * \return PanNDE::CArray<T> Structure containing pointer to data and its length
      * \warning The returned structure does not own the data - ensure proper lifetime management
      */
      const PanNDE::CArray<T> getCArray()override{
        PanNDE::CArray<T> carray;
        carray.data=vec.data();
        carray.length=vec.size();
        return carray;
      }

    private:
      //! Internal vector for storing array elements
      std::vector<T> vec;
  };

  /*! \class HostDataBundle<T> HostArray.hpp "modules/HostData/include/HostArray.hpp"
  *
  * Provides a container for grouping named scalars and arrays together.
  * This class implements the PanNDE::DataBundle interface, allowing users to
  * store related data elements with string identifiers for easy retrieval.
  *
  */
  template<typename T>
  class HostDataBundle : public PanNDE::DataBundle<T> {
    public:
      /*!
      * Retrieves a scalar value by its name.
      * \param keyname std::string The name of the scalar to retrieve
      * \return T The scalar value
      * \throw std::out_of_range If a scalar with the given name doesn't exist
      */
      T scalar(std::string keyname)override{return scalars.at(keyname);};
      
      /*!
      * Retrieves an array by its name.
      * \param keyname std::string The name of the array to retrieve
      * \return std::shared_ptr<PanNDE::Array<T>> Shared pointer to the requested array
      * \throw std::out_of_range If an array with the given name doesn't exist
      */
      std::shared_ptr<PanNDE::Array<T>> array(std::string keyname)override{
        return arrays.at(keyname);
      };
      
      /*!
      * Gets the name of a scalar by its index.
      * \param idx int Index of the scalar name to retrieve
      * \return std::string The name of the scalar at the specified index
      * \throw std::out_of_range If the index is out of bounds
      */
      std::string scalarName(int idx)override{return scalarnames.at(idx);};
      
      /*!
      * Gets the name of an array by its index.
      * \param idx int Index of the array name to retrieve
      * \return std::string The name of the array at the specified index
      * \throw std::out_of_range If the index is out of bounds
      */
      std::string arrayName(int idx)override{return arraynames.at(idx);};
      
      /*!
      * Returns the number of scalars in the bundle.
      * \return int The number of scalar values stored in this bundle
      */
      int scalarCount()override{return scalarnames.size();};
      
      /*!
      * Returns the number of arrays in the bundle.
      * \return int The number of arrays stored in this bundle
      */
      int arrayCount()override{return arraynames.size();};
      
      /*!
      * Adds a new scalar to the bundle with the given name.
      * \param keyname std::string Name to associate with the scalar
      * \param value T The scalar value to store
      * \throw std::runtime_error If a scalar with the given name already exists
      */
      void emplaceScalar(std::string keyname,T value)override{
        auto result=scalars.emplace(keyname,value);
        addScalarKey(result.second,keyname);
      };
      
      /*!
      * Adds a new array to the bundle with the given name.
      * \param keyname std::string Name to associate with the array
      * \param array std::shared_ptr<PanNDE::Array<T>> The array to store
      * \throw std::runtime_error If an array with the given name already exists
      */
      void emplaceArray(std::string keyname,std::shared_ptr<PanNDE::Array<T>> array)override{
        auto result=arrays.emplace(keyname,array);
        addArrayKey(result.second,keyname);
      };
      
    private:
      /*!
      * Helper method to register a scalar name after insertion attempt.
      * \param inserted bool Whether the insertion was successful
      * \param keyname std::string The name that was used for insertion
      * \throw std::runtime_error If the scalar name is already in use
      */
      void addScalarKey(bool inserted,std::string keyname){
        if(inserted){
          scalarnames.push_back(keyname);
        }else{throw std::runtime_error("scalar " + keyname + " already in use.");};
      };
      
      /*!
      * Helper method to register an array name after insertion attempt.
      * \param inserted bool Whether the insertion was successful
      * \param keyname std::string The name that was used for insertion
      * \throw std::runtime_error If the array name is already in use
      */
      void addArrayKey(bool inserted,std::string keyname){
        if(inserted){
          arraynames.push_back(keyname);
        }else{throw std::runtime_error("array " + keyname + " already in use.");};
      };
      
      //! Stores the names of arrays in insertion order
      std::vector<std::string> arraynames;
      
      //! Stores the names of scalars in insertion order
      std::vector<std::string> scalarnames;
      
      //! Maps array names to their corresponding array objects
      std::map<std::string,std::shared_ptr<PanNDE::Array<T>>> arrays;
      
      //! Maps scalar names to their corresponding values
      std::map<std::string,T> scalars;
  };

  /*! \class HostArrayFactory<T> HostArray.hpp "modules/HostData/include/HostArray.hpp"
  * 
  * Factory class for creating host-based array and data bundle objects.
  * Implements the PanNDE::ArrayFactory interface to provide a consistent
  * mechanism for creating array objects appropriate for the host environment.
  *
  */
  template<typename T>
  class HostArrayFactory : public PanNDE::ArrayFactory<T> {
    public:
      /*!
      * Creates a shared pointer to a new HostArrayFactory instance.
      * \return std::shared_ptr<HostData::HostArrayFactory<T>> Shared pointer to the new factory
      */
      static
      std::shared_ptr<HostData::HostArrayFactory<T>> makeShared(){
        auto factory=std::make_shared<HostData::HostArrayFactory<T>>(HostData::HostArrayFactory<T>());
        return std::move(factory);
      };
      
      /*!
      * Creates a new empty managed array.
      * \return std::shared_ptr<PanNDE::Array<T>> Shared pointer to the newly created array
      */
      std::shared_ptr<PanNDE::Array<T>> makeManagedArray()override{
        std::shared_ptr<PanNDE::Array<T>> array;
        array=std::make_shared<HostData::HostArray<T>>(HostData::HostArray<T>());
        return std::move(array);
      };
      
      /*!
      * Creates a new raw array object (caller takes ownership).
      * \return PanNDE::Array<T>* Pointer to the newly created array
      * \warning Use with caution - caller is responsible for memory management
      */
      PanNDE::Array<T>* newArray()override{
        auto array=new HostData::HostArray<T>();
        return (PanNDE::Array<T>*)array;
      };
      
      /*!
      * Deletes an array created by newArray().
      * \param array PanNDE::Array<T>* Pointer to the array to delete
      * \warning Only use with arrays created by this factory's newArray() method
      */
      void deleteArray(PanNDE::Array<T>* array)override{delete array;};

      /*!
      * Creates a new empty managed data bundle.
      * \return std::shared_ptr<PanNDE::DataBundle<T>> Shared pointer to the newly created bundle
      */
      std::shared_ptr<PanNDE::DataBundle<T>> makeManagedDataBundle()override{
        auto bundle=std::make_shared<HostData::HostDataBundle<T>>(HostData::HostDataBundle<T>());
        return std::move(bundle);
      };
  };
};