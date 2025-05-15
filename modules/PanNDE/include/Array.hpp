/*! \headerfile Array.hpp "modules/PanNDE/include/Array.hpp"
* "Array.hpp" defines the core data container interfaces for PanNDE.
* It provides platform-independent array abstractions that can work across
* various hardware architectures (CPU, GPU, etc.) without relying on
* standard library containers. This design enables implementation flexibility
* while maintaining consistent interfaces for simulation components.
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
  /*! \struct CArray
  * 
  * A lightweight C-compatible array structure for low-level data transfer.
  * Provides a simple representation of an array with a pointer and length
  * that can be used for interoperability with external libraries and systems.
  *
  * \tparam T The data type of the array elements
  */
  template<typename T>
  struct CArray {
    T* data=nullptr;       //!< Pointer to the raw data array
    int64_t length=0;      //!< Number of elements in the array
  };

  /*! \class Array<T> Array.hpp "modules/PanNDE/include/Array.hpp"
  *
  * Defines a platform-independent interface for fixed-size arrays.
  * 
  * This interface provides a common set of operations for array manipulation
  * without depending on specific implementations like STL containers. This
  * abstraction enables implementations for different hardware architectures
  * including CPUs and GPUs where STL might be unavailable or inefficient.
  *
  * Implementations of this interface should focus on:
  * - Efficient memory management for the target architecture
  * - Safe element access with bounds checking
  * - Appropriate resize behavior for the intended use case
  *
  * \tparam T The data type of the array elements
  */
  template<typename T>
  class Array {
    public:
      /*! 
      * Access an array element by reference.
      * 
      * Implementations should provide bounds checking and throw exceptions
      * for out-of-range accesses.
      * 
      * \param idx The zero-based index of the element to access
      * \return T& Reference to the requested element
      * \throws Exception if idx is out of range
      */
      virtual T& at(int64_t idx) =0;
      
      /*! 
      * Get the number of elements in the array.
      * 
      * \return int64_t The current size of the array
      */
      virtual int64_t size() =0;
      
      /*! 
      * Get direct access to the underlying array data.
      * 
      * This method provides raw pointer access for performance-critical
      * operations or interfacing with external libraries. Use with caution
      * as it bypasses safety features like bounds checking.
      * 
      * \return T* Pointer to the first element of the array
      */
      virtual T* data() =0;

      /*! 
      * Resize the array to a new size.
      * 
      * Implementations may choose whether to preserve existing elements
      * when resizing. The behavior should be documented by each implementation.
      * 
      * \param newsize The new size for the array
      */
      virtual void resize(int64_t newsize) =0;
      
      /*! 
      * Copy data from another array using a raw pointer.
      * 
      * \param other Pointer to the source array to copy from
      */
      virtual void copy(PanNDE::Array<T>* other) =0;
      
      /*! 
      * Copy data from another array using a shared pointer.
      * 
      * \param other Shared pointer to the source array to copy from
      */
      virtual void copy(std::shared_ptr<PanNDE::Array<T>> other) =0;

      /*! 
      * Copy data from a CArray structure.
      * 
      * This method is primarily intended for interprocess communication
      * or integration with external systems. CArray provides no bounds
      * protection, so this method should be used with caution.
      * 
      * \param other Reference to a CArray to copy from
      */
      virtual void copy(CArray<T>& other) =0;
      
      /*! 
      * Export array data as a CArray structure.
      * 
      * Provides a lightweight representation of the array for low-level
      * operations or data exchange. The returned CArray references the
      * original data and does not make a copy.
      * 
      * \return CArray<T> A structure containing a pointer to the array data and its size
      */
      virtual const CArray<T> getCArray() =0;

      //! Virtual destructor to ensure proper cleanup in derived classes
      virtual ~Array(){};
  };

  /*! \class DataBundle<T> Array.hpp "modules/PanNDE/include/Array.hpp"
  *
  * A container for named scalars and arrays of a common data type.
  * 
  * DataBundle provides a way to group related data items and access them by name.
  * It serves as a flexible container for both scalar values and arrays, allowing
  * simulation components to exchange structured data collections.
  *
  * Typical uses include:
  * - Collecting input parameters for simulations
  * - Bundling multiple result arrays from an analysis
  * - Passing metadata along with primary data
  *
  * \tparam T The data type for all values and arrays in the bundle
  */
  template<typename T>
  class DataBundle {
    public:
      /*! 
      * Get a scalar value by its name.
      * 
      * \param keyname The name identifier of the scalar to retrieve
      * \return T The scalar value
      * \throws Exception if the name does not exist in the bundle
      */
      virtual T scalar(std::string keyname) =0;
      
      /*! 
      * Get an array by its name.
      * 
      * \param keyname The name identifier of the array to retrieve
      * \return std::shared_ptr<PanNDE::Array<T>> Shared pointer to the requested array
      * \throws Exception if the name does not exist in the bundle
      */
      virtual std::shared_ptr<PanNDE::Array<T>> array(std::string keyname) =0;
      
      /*! 
      * Get the name of a scalar by its index.
      * 
      * \param idx The index of the scalar name to retrieve
      * \return std::string The name of the scalar at the given index
      * \throws Exception if the index is out of range
      */
      virtual std::string scalarName(int idx) =0;
      
      /*! 
      * Get the name of an array by its index.
      * 
      * \param idx The index of the array name to retrieve
      * \return std::string The name of the array at the given index
      * \throws Exception if the index is out of range
      */
      virtual std::string arrayName(int idx) =0;
      
      /*! 
      * Get the number of scalars in the bundle.
      * 
      * \return int The number of scalar values
      */
      virtual int scalarCount() =0;
      
      /*! 
      * Get the number of arrays in the bundle.
      * 
      * \return int The number of arrays
      */
      virtual int arrayCount() =0;
      
      /*! 
      * Add a scalar value to the bundle with a specified name.
      * 
      * \param keyname The name to associate with the scalar
      * \param value The scalar value to add
      * \throws Exception if the name already exists in the bundle
      */
      virtual void emplaceScalar(std::string keyname, T value) =0;
      
      /*! 
      * Add an array to the bundle with a specified name.
      * 
      * \param keyname The name to associate with the array
      * \param array Shared pointer to the array to add
      * \throws Exception if the name already exists in the bundle
      */
      virtual void emplaceArray(std::string keyname,
                                std::shared_ptr<PanNDE::Array<T>> array) =0;
  };

  /*! \class ArrayFactory<T> Array.hpp "modules/PanNDE/include/Array.hpp"
  *
  * Factory interface for creating Array and DataBundle instances.
  * 
  * This interface provides a standardized way to instantiate platform-specific
  * implementations of the Array and DataBundle interfaces. It encapsulates the
  * creation logic, allowing simulation components to work with these data
  * structures without knowing the specific implementation details.
  *
  * Implementing this factory for different platforms (CPU, GPU, etc.) allows
  * the same simulation code to run across different hardware architectures.
  *
  * \tparam T The data type for the arrays and data bundles to be created
  */
  template<typename T>
  class ArrayFactory {
    public:
      /*! 
      * Create a new empty array managed by a shared pointer.
      * 
      * This is the preferred method for creating arrays as it provides
      * automatic memory management.
      * 
      * \return std::shared_ptr<PanNDE::Array<T>> A shared pointer to a newly created array
      */
      virtual std::shared_ptr<PanNDE::Array<T>> makeManagedArray() =0;
      
      /*! 
      * Create a new empty array with manual memory management.
      * 
      * This method requires manual deletion using deleteArray().
      * It's provided primarily for compatibility with C-style interfaces.
      * 
      * \return PanNDE::Array<T>* A raw pointer to a newly created array
      */
      virtual PanNDE::Array<T>* newArray() =0;
      
      /*! 
      * Delete an array created with newArray().
      * 
      * \param array Pointer to the array to delete
      */
      virtual void deleteArray(PanNDE::Array<T>* array) =0;
      
      /*! 
      * Create a new empty data bundle managed by a shared pointer.
      * 
      * \return std::shared_ptr<PanNDE::DataBundle<T>> A shared pointer to a newly created data bundle
      */
      virtual std::shared_ptr<PanNDE::DataBundle<T>> makeManagedDataBundle() =0;
  };
};