/*! \headerfile Field.hpp "modules/PanNDE/include/Field.hpp"
* "Field.hpp" defines interfaces for simulation data fields that are associated with
* computational meshes. Fields represent physical quantities or solution variables
* distributed across a mesh, either at node points or cell centers.
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

#include "Mesh.hpp"

namespace PanNDE{
  /*! \class Field Field.hpp "modules/PanNDE/include/Field.hpp"
  *
  * Represents a physical quantity or solution variable distributed across a mesh.
  *
  * Fields associate scalar values with either mesh nodes or cells, providing a way
  * to represent continuous physical quantities in the discrete computational domain.
  * Common field examples include displacement, stress, temperature, and material
  * properties.
  *
  * Fields are strongly tied to a specific mesh and provide methods for accessing, 
  * manipulating, and interpolating values between different representation types
  * (e.g., node-based vs. cell-based).
  *
  */
  class Field{
    public:
      /*! 
      * Enumeration defining where field data is stored in relation to the mesh.
      */
      enum FieldType{
        NODE=0,     //!< Values stored at mesh nodes (vertices)
        CELL=1,     //!< Values stored at cell centers
        CONSTANT=2  //!< Single value applied to the entire field (treated as cell-centered internally)
      };

      /*!
      * Gets the storage type of this field.
      *
      * The field type determines how the data values are associated with the mesh
      * (at nodes, cells, or constant across the domain) and affects interpolation
      * behavior between different field types.
      *
      * \return FieldType The storage type (NODE, CELL, or CONSTANT)
      */
      virtual FieldType type() =0;

      /*!
      * Gets the mesh to which this field is attached.
      *
      * Every field is associated with a specific mesh that defines the spatial
      * domain on which the field values are distributed.
      *
      * \return std::shared_ptr<PanNDE::Mesh> Shared pointer to the associated mesh
      */
      virtual std::shared_ptr<PanNDE::Mesh> mesh() =0;

      /*!
      * Gets direct access to the underlying data array.
      *
      * This method provides low-level access to the raw field data and should be
      * used with caution, primarily for performance-critical operations or when
      * interfacing with external libraries that require direct memory access.
      *
      * \return double* Pointer to the first element of the data array
      */
      virtual double* data() =0;

      /*!
      * Gets a reference to a field value at the specified index.
      *
      * The index corresponds to either a node or cell index, depending on the
      * field's type. This provides a way to modify values in-place.
      *
      * \param index int Location index to access (node index for NODE fields, cell index for CELL fields)
      * \return double& Reference to the value at the specified index
      */
      virtual double& at(int index) =0;

      /*!
      * Gets the field value at a specified cell.
      *
      * For CELL and CONSTANT fields, this directly retrieves the value.
      * For NODE fields, this may involve interpolation from nodal values.
      *
      * \param cidx int Cell index
      * \return double The field value at the specified cell
      */
      virtual double atCell(int cidx) =0;

      /*!
      * Gets the field value at a specified node.
      *
      * For NODE fields, this directly retrieves the value.
      * For CELL and CONSTANT fields, this may involve interpolation from cell values.
      *
      * \param nidx int Node index
      * \return double The field value at the specified node
      */
      virtual double atNode(int nidx) =0;

      /*!
      * Gets the number of values in the field.
      *
      * This returns the total count of values stored in the field, which corresponds
      * to either the number of nodes or cells in the mesh, depending on the field type.
      * For CONSTANT fields, this may return 1.
      *
      * \return int32_t Number of field values
      */
      virtual int32_t size() =0;

      /*!
      * Maps data from another field into this field.
      *
      * If the source and destination fields are of the same type (e.g., both NODE),
      * this performs a direct copy. If they differ (e.g., NODE to CELL), this performs
      * an appropriate interpolation to transform the data.
      *
      * Both fields should be associated with the same mesh.
      *
      * \param other std::shared_ptr<PanNDE::Field> Source field to map from
      */
      virtual void mapFrom(std::shared_ptr<PanNDE::Field> other) =0;

      /*!
      * Maps data from another field into this field using a raw pointer.
      *
      * This is an alternative to the shared pointer version, with the same behavior.
      * Not recommended for general use as it doesn't leverage memory safety features.
      *
      * \param other PanNDE::Field* Source field to map from
      */
      virtual void mapFrom(PanNDE::Field* other) =0;

      //! Virtual destructor to ensure proper cleanup in derived classes
      virtual ~Field(){};
  };

  /*! \class FieldBundle Field.hpp "modules/PanNDE/include/Field.hpp"
  *
  * A collection of related fields associated with a common mesh.
  *
  * FieldBundle provides a way to group multiple fields that represent different
  * aspects of a physical system or solution. All fields in a bundle share the
  * same underlying mesh, ensuring consistency in the spatial domain.
  *
  * Typical uses include:
  * - Grouping component fields of a vector quantity (e.g., displacement x,y,z)
  * - Collecting various solution variables (e.g., stress, strain, temperature)
  * - Managing a complete solution state for a simulation
  *
  */
  class FieldBundle {
    public:
      /*!
      * Gets the shared mesh that all fields in this bundle are associated with.
      *
      * \return std::shared_ptr<PanNDE::Mesh>& Reference to the shared mesh pointer
      */
      virtual std::shared_ptr<PanNDE::Mesh>& mesh() =0;

      /*!
      * Gets a field from the bundle by its name.
      *
      * \param keyname std::string Name of the field to retrieve
      * \return std::shared_ptr<PanNDE::Field> The requested field
      * \throws Exception if the field name doesn't exist in the bundle
      */
      virtual std::shared_ptr<PanNDE::Field> field(std::string keyname) =0;

      /*!
      * Gets the name of a field by its index in the bundle.
      *
      * \param idx int Index of the field
      * \return std::string Name of the field at the specified index
      * \throws Exception if the index is out of range
      */
      virtual std::string fieldName(int idx) =0;

      /*!
      * Gets the total number of fields in the bundle.
      *
      * \return int Number of fields
      */
      virtual int fieldCount() =0;

      /*!
      * Adds an existing field to the bundle with the specified name.
      *
      * The field must be associated with the same mesh as other fields in the bundle.
      *
      * \param keyname std::string Name to assign to the field
      * \param field std::shared_ptr<PanNDE::Field> The field to add
      * \throws Exception if a field with the same name already exists or if the mesh doesn't match
      */
      virtual void emplaceField(std::string keyname,std::shared_ptr<PanNDE::Field> field) =0;

      /*!
      * Creates and adds a new field to the bundle with the specified name and type.
      *
      * The new field will be associated with the bundle's mesh.
      *
      * \param keyname std::string Name to assign to the field
      * \param type PanNDE::Field::FieldType Type of field to create (NODE, CELL, or CONSTANT)
      * \throws Exception if a field with the same name already exists
      */
      virtual void emplaceField(std::string keyname,PanNDE::Field::FieldType type) =0;

      //! Virtual destructor to ensure proper cleanup in derived classes
      virtual ~FieldBundle(){};
  };

  /*! \class FieldFactory Field.hpp "modules/PanNDE/include/Field.hpp"
  *
  * Factory interface for creating fields and field collections.
  *
  * This interface provides methods to create individual fields and field bundles
  * with different memory management strategies. Implementations should create fields
  * that are properly associated with the specified mesh and initialized according
  * to the requested field type.
  *
  */
  class FieldFactory{
    public:
      /*!
      * Creates a field with shared pointer management.
      *
      * This is the preferred method for creating fields as it provides automatic
      * memory management.
      *
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to associate with the field
      * \param tp PanNDE::Field::FieldType The type of field to create
      * \return std::shared_ptr<PanNDE::Field> A shared pointer to the newly created field
      */
      virtual std::shared_ptr<PanNDE::Field> makeManagedField(std::shared_ptr<PanNDE::Mesh> mesh,
                                                              PanNDE::Field::FieldType tp) =0;

      /*!
      * Creates a field with shared pointer management using a raw mesh pointer.
      *
      * \param mesh PanNDE::Mesh* The mesh to associate with the field
      * \param tp PanNDE::Field::FieldType The type of field to create
      * \return std::shared_ptr<PanNDE::Field> A shared pointer to the newly created field
      */
      virtual std::shared_ptr<PanNDE::Field> makeManagedField(PanNDE::Mesh* mesh,
                                                              PanNDE::Field::FieldType tp) =0;

      /*!
      * Creates a field with manual memory management.
      *
      * This method is not recommended for general use and requires manual deletion
      * using deleteField(). It's provided primarily for compatibility with C-style interfaces.
      *
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to associate with the field
      * \param tp PanNDE::Field::FieldType The type of field to create
      * \return PanNDE::Field* A raw pointer to the newly created field
      */
      virtual PanNDE::Field* newField(std::shared_ptr<PanNDE::Mesh> mesh,
                                      PanNDE::Field::FieldType tp) =0;

      /*!
      * Creates a field with manual memory management using a raw mesh pointer.
      *
      * This method is not recommended for general use and requires manual deletion
      * using deleteField(). It's provided primarily for compatibility with C-style interfaces.
      *
      * \param mesh PanNDE::Mesh* The mesh to associate with the field
      * \param tp PanNDE::Field::FieldType The type of field to create
      * \return PanNDE::Field* A raw pointer to the newly created field
      */
      virtual PanNDE::Field* newField(PanNDE::Mesh* mesh,
                                      PanNDE::Field::FieldType tp) =0;

      /*!
      * Deletes a field that was created using newField().
      *
      * This method must be called for any field created with newField() to prevent memory leaks.
      *
      * \param field PanNDE::Field* Pointer to the field to delete
      */
      virtual void deleteField(PanNDE::Field* field) =0;
      
      /*!
      * Creates an array to store multiple fields.
      *
      * This is useful when working with collections of related fields that
      * need to be processed together but don't necessarily share the same mesh.
      *
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> An array of field pointers
      */
      virtual std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> makeManagedFieldArray() =0;

      /*!
      * Creates an empty field bundle.
      *
      * The bundle can then be populated with fields using the emplaceField methods.
      *
      * \return std::shared_ptr<PanNDE::FieldBundle> A shared pointer to the newly created field bundle
      */
      virtual std::shared_ptr<PanNDE::FieldBundle> makeEmptyManagedFieldBundle() =0;
  };
};