/*! \headerfile HostArray.hpp "modules/HostData/include/HostArray.hpp"
* "HostArray.hpp" contains the class implementation encapsulating 
* a fixed size array. This allows a "in-family" transfer of arrays
* without requiring the use of a library, including the STL. The 
* rationale being that extention to GPU or other device code usually 
* prohibits the STL, and therefore avoidance of using std::vector<T> or 
* std::array<T,int> in interface definitions becomes a requirement to 
* ensure interface extensibility.
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
  * Implements the expected characteristics of a fixed-size array class. 
  * This class allows a "in-family" transfer of arrays without requiring 
  * the use of a library, including the STL. The rationale being that extention 
  * to GPU or other device code usually prohibits the STL, and therefore 
  * avoidance of using std::vector<T> or std::array<T,int> in interface 
  * definitions becomes a requirement to ensure interface extensibility.
  *
  */
  template<typename T>
  class HostArray : public PanNDE::Array<T> {
    public:
      /*! 
      * Get array value by reference.
      * \param idx array index to access
      */
      T& at(int64_t idx)override{return vec.at(idx);};
      /*! 
      * Get array size.
      */
      int64_t size()override{return vec.size();};
      /*! 
      * Resize array. Array contents not guaranteed to be preserved
      * \param newsize new size of array.
      */
      void resize(int64_t newsize)override{return vec.resize(newsize);};
      /*! 
      * Get C-pointer to head of array. Only use when low-level access is required.
      */
      T* data()override{return vec.data();};

      /*! 
      * Copy other array to the local array
      * \param other C-pointer to source array for copy.
      */
      void copy(PanNDE::Array<T>* other)override{
        this->resize(other->size());
        for(auto k=0;k<(this->size());k++){this->at(k)=other->at(k);};
      };
      /*! 
      * Copy other array to the local array
      * \param other shared pointer to source array for copy.
      */
      void copy(std::shared_ptr<PanNDE::Array<T>> other)override{
        copy(other.get());
      };

      /*! 
      * Copy other array to the local array. This is not recommended for general use, 
      * CArray has no overrun protections. This is primarily for inter-process or device 
      * communications where low-level access is required.
      * \param other reference to Carray wrapped struct for copy.
      */
      void copy(PanNDE::CArray<T>& other)override{
        this->resize(other.length);
        for(auto k=0;k<(this->size());k++){this->at(k)=other.data[k];};
      };
      /*! 
      * Get array data as CArray struct. This is not recommended for general use, 
      * CArray has no overrun protections. This is primarily for inter-process or device 
      * communications where low-level access is required.
      */
      const PanNDE::CArray<T> getCArray()override{
        PanNDE::CArray<T> carray;
        carray.data=vec.data();
        carray.length=vec.size();
        return carray;
      }

    private:
      std::vector<T> vec;
  };

  /*! \class HostDataBundle<T> HostArray.hpp "modules/HostData/include/HostArray.hpp"
  *
  * Implements a class which bundles scalars and arrays by name. 
  * This can be used to group loosely associated data 
  *
  */
  template<typename T>
  class HostDataBundle : public PanNDE::DataBundle<T> {
    public:
      /*! 
      * get scalar value by name
      */
      T scalar(std::string keyname)override{return scalars.at(keyname);};
      /*! 
      * get array by name
      */
      std::shared_ptr<PanNDE::Array<T>> array(std::string keyname)override{
        return arrays.at(keyname);
      };
      /*! 
      * get scalar name by index
      */
      std::string scalarName(int idx)override{return scalarnames.at(idx);};
      /*! 
      * get array name by index
      */
      std::string arrayName(int idx)override{return arraynames.at(idx);};
      /*! 
      * get number of scalars
      */
      int scalarCount()override{return scalarnames.size();};
      /*! 
      * get number of arrays
      */
      int arrayCount()override{return arraynames.size();};
      /*! 
      * add scalar to bundle
      */
      void emplaceScalar(std::string keyname,T value)override{
        auto result=scalars.emplace(keyname,value);
        addScalarKey(result.second,keyname);
      };
      /*! 
      * add array to bundle
      */
      void emplaceArray(std::string keyname,std::shared_ptr<PanNDE::Array<T>> array)override{
        auto result=arrays.emplace(keyname,array);
        addArrayKey(result.second,keyname);
      };
      
    private:
      void addScalarKey(bool inserted,std::string keyname){
        if(inserted){
          scalarnames.push_back(keyname);
        }else{throw std::runtime_error("scalar " + keyname + " already in use.");};
      };
      void addArrayKey(bool inserted,std::string keyname){
        if(inserted){
          arraynames.push_back(keyname);
        }else{throw std::runtime_error("array " + keyname + " already in use.");};
      };
      std::vector<std::string> arraynames;
      std::vector<std::string> scalarnames;
      std::map<std::string,std::shared_ptr<PanNDE::Array<T>>> arrays;
      std::map<std::string,T> scalars;
  };

  /*! \class HostArrayFactory<T> HostArray.hpp "modules/HostData/include/HostArray.hpp"
  * Implements a factory class to create the HostArray<T> class
  *
  */
  template<typename T>
  class HostArrayFactory : public PanNDE::ArrayFactory<T> {
    public:
      /*! 
      * create an empty shared array
      */
      std::shared_ptr<PanNDE::Array<T>> makeManagedArray()override{
        std::shared_ptr<PanNDE::Array<T>> array;
        array=std::make_shared<HostData::HostArray<T>>(HostData::HostArray<T>());
        return std::move(array);
      };
      /*! 
      * create an empty array. Not recommended, but included for the use case
      */
      PanNDE::Array<T>* newArray()override{
        auto array=new HostData::HostArray<T>();
        return (PanNDE::Array<T>*)array;
      };
      /*! 
      * delete an array created with newArray(). Not recommended, but included for the use case
      */
      void deleteArray(PanNDE::Array<T>* array)override{delete array;};

      /*! 
      * create an empty shared data bundle
      */
      std::shared_ptr<PanNDE::DataBundle<T>> makeManagedDataBundle()override{
        auto bundle=std::make_shared<HostData::HostDataBundle<T>>(HostData::HostDataBundle<T>());
        return std::move(bundle);
      };
  };
};