/*! \headerfile Array.hpp "modules/PanNDE/include/Array.hpp"
* "Array.hpp" contains the class definition encapsulating 
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

namespace PanNDE {
  template<typename T>
  struct CArray {
    T* data=nullptr;
    int64_t length=0;
  };

  /*! \class Array<T> Array.hpp "modules/PanNDE/include/Array.hpp"
  *
  * Defines the expected characteristics of a fixed-size array class. 
  * This class allows a "in-family" transfer of arrays without requiring 
  * the use of a library, including the STL. The rationale being that extention 
  * to GPU or other device code usually prohibits the STL, and therefore 
  * avoidance of using std::vector<T> or std::array<T,int> in interface 
  * definitions becomes a requirement to ensure interface extensibility.
  *
  */
  template<typename T>
  class Array {
    public:
      /*! 
      * Get array value by reference.
      * \param idx array index to access
      */
      virtual T& at(int64_t idx) =0;
      /*! 
      * Get array size.
      */
      virtual int64_t size() =0;
      /*! 
      * Get C-pointer to head of array. Only use when low-level access is required.
      */
      virtual T* data() =0;

      /*! 
      * Resize array. Array contents not guaranteed to be preserved
      * \param newsize new size of array.
      */
      virtual void resize(int64_t newsize) =0;
      
      /*! 
      * Copy other array to the local array
      * \param other C-pointer to source array for copy.
      */
      virtual void copy(PanNDE::Array<T>* other) =0;
      /*! 
      * Copy other array to the local array
      * \param other shared pointer to source array for copy.
      */
      virtual void copy(std::shared_ptr<PanNDE::Array<T>> other) =0;

      /*! 
      * Copy other array to the local array. This is not recommended for general use, 
      * CArray has no overrun protections. This is primarily for inter-process or device 
      * communications where low-level access is required.
      * \param other reference to Carray wrapped struct for copy.
      */
      virtual void copy(CArray<T>& other) =0;
      /*! 
      * Get array data as CArray struct. This is not recommended for general use, 
      * CArray has no overrun protections. This is primarily for inter-process or device 
      * communications where low-level access is required.
      */
      virtual const CArray<T> getCArray() =0;
  };

  /*! \class DataBundle<T> Array.hpp "modules/PanNDE/include/Array.hpp"
  *
  * Defines a class which bundles scalars and arrays by name. 
  * This can be used to group loosely associated data 
  *
  */
  template<typename T>
  class DataBundle {
    public:
      /*! 
      * get scalar value by name
      */
      virtual T scalar(std::string kyename) =0;
      /*! 
      * get array by name
      */
      virtual std::shared_ptr<PanNDE::Array<T>> array(std::string keyname) =0;
      /*! 
      * get scalar name by index
      */
      virtual std::string scalarName(int idx) =0;
      /*! 
      * get array name by index
      */
      virtual std::string arrayName(int idx) =0;
      /*! 
      * get number of scalars
      */
      virtual int scalarCount() =0;
      /*! 
      * get number of arrays
      */
      virtual int arrayCount() =0;
      /*! 
      * add scalar to bundle
      */
      virtual void emplaceScalar(std::string keyname,T value) =0;
      /*! 
      * add array to bundle
      */
      virtual void emplaceArray(std::string keyname,
                                std::shared_ptr<PanNDE::Array<T>> array) =0;
  };

  /*! \class ArrayFactory<T> Array.hpp "modules/PanNDE/include/Array.hpp"
  *
  * Defines a factory class to create the Array<T> class
  *
  */
  template<typename T>
  class ArrayFactory {
    public:
      /*! 
      * create an empty shared array
      */
      virtual std::shared_ptr<PanNDE::Array<T>> makeManagedArray() =0;
      /*! 
      * create an empty array. Not recommended, but included for the use case
      */
      virtual PanNDE::Array<T>* newArray() =0;
      /*! 
      * delete an array created with newArray(). Not recommended, but included for the use case
      */
      virtual void deleteArray(PanNDE::Array<T>* array) =0;
      /*! 
      * create an empty shared data bundle
      */
      virtual std::shared_ptr<PanNDE::DataBundle<T>> makeManagedDataBundle() =0;
  };
};