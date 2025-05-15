/*! \headerfile HostField.hpp "modules/HostData/include/HostField.hpp"
* "HostField.hpp" contains class implementations for field data attached to meshes.
* It provides host-based data storage for scalar fields defined on nodes or cells,
* field bundles for organizing multiple related fields, and factory methods for
* field creation.
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

#include <vector>
#include <map>

#include "Array.hpp"
#include "Field.hpp"
#include "Mesh.hpp"

#include "HostArray.hpp"

namespace HostData {
  /*! \class HostField HostField.hpp "modules/HostData/include/HostField.hpp"
  * 
  * Implements a scalar field attached to a mesh.
  * This class stores numerical data associated with a mesh, either at nodes,
  * cells, or as a constant value. It provides methods to access values
  * directly or to interpolate between node and cell values as needed.
  *
  */
  class HostField : public PanNDE::Field {
    public:
      /*!
      * Enumeration defining the possible field data locations.
      */
      enum FieldType{
        NODE=0,      //!< Field data stored at mesh nodes (vertices)
        CELL=1,      //!< Field data stored at mesh cells (elements)
        CONSTANT=2   //!< Single value for the entire field (treated as cell-centered when interpolating)
      };
      
      /*!
      * Gets the storage type of this field.
      * \return PanNDE::Field::FieldType The type of field (NODE, CELL, or CONSTANT)
      */
      PanNDE::Field::FieldType type()override{return mytype;};
      
      /*!
      * Gets the mesh to which this field is attached.
      * \return std::shared_ptr<PanNDE::Mesh> The mesh associated with this field
      */
      std::shared_ptr<PanNDE::Mesh> mesh()override{return themesh;};
      
      /*!
      * Gets a raw pointer to the underlying field data.
      * \return double* Pointer to the internal data array
      * \warning Use with caution, only when low-level direct access is necessary
      */
      double* data()override{return field_data.data();};
      
      /*!
      * Gets a reference to the field value at the specified index.
      * For CONSTANT fields, always returns a reference to the single value.
      * \param index int Index corresponding to a node or cell ID (based on field type)
      * \return double& Reference to the field value
      * \throw std::out_of_range If index is invalid
      */
      double& at(int index)override{
        return ((PanNDE::Field::CONSTANT==mytype)?field_data.at(0):field_data.at(index));
      };

      /*!
      * Gets the field value at the specified cell.
      * For NODE fields, averages values from the cell's nodes.
      * For CELL or CONSTANT fields, returns the direct value.
      * \param cidx int Cell index
      * \return double The field value at the specified cell
      */
      double atCell(int cidx)override{
        if(PanNDE::Field::NODE==mytype){
          auto box=themesh->cell(cidx);
          double value=0.0;
          for(int k=0;k<box->size();k++){
            value+=(-1==box->at(k))?0.:at(box->at(k));
            //count+=(-1==box->at(k))?0:1;
          };
          value=value/(double(box->size()));
          return value;
        };
        return at(cidx);
      };
      
      /*!
      * Gets the field value at the specified node.
      * For CELL fields, averages values from cells connected to the node.
      * For NODE or CONSTANT fields, returns the direct value.
      * \param nidx int Node index
      * \return double The field value at the specified node
      */
      double atNode(int nidx)override{
        if(PanNDE::Field::NODE==mytype){
          return at(nidx);
        };
        auto box=themesh->connectedCells(nidx);
        double value=0.0;//int count=0;
        for(int k=0;k<box->size();k++){
          value+=(-1==box->at(k))?0.:at(box->at(k));
          //count+=(-1==box->at(k))?0:1;
        };
        value=value/(double(box->size()));
        return value;
      };

      /*!
      * Gets the number of values in the field data array.
      * \return int32_t Number of stored values (equal to node count, cell count, or 1)
      */
      int32_t size()override{return field_data.size();};

      /*!
      * Maps data from another field into this field.
      * If field types match, performs a direct copy.
      * If field types differ, performs appropriate interpolation.
      * \param other std::shared_ptr<PanNDE::Field> The source field
      */
      void mapFrom(std::shared_ptr<PanNDE::Field> other)override{
        mapFrom(other.get());
      };
      
      /*!
      * Maps data from another field into this field.
      * If field types match, performs a direct copy.
      * If field types differ, performs appropriate interpolation.
      * \param other PanNDE::Field* Raw pointer to the source field
      */
      void mapFrom(PanNDE::Field* other)override{
        themesh=other->mesh();
        if(PanNDE::Field::CONSTANT==mytype){setUp(other->type());};
        for(auto k=0;k<size();k++){
          at(k)=(PanNDE::Field::NODE==mytype)?other->atNode(k):other->atCell(k);
        };
      };

      /*!
      * Constructor for creating a field attached to a mesh.
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to attach this field to
      * \param tp PanNDE::Field::FieldType The storage type of this field
      * \throw std::runtime_error If field type is not supported
      */
      HostField(std::shared_ptr<PanNDE::Mesh> mesh,PanNDE::Field::FieldType tp){
        themesh=mesh;
        setUp(tp);
      };
      
      /*!
      * Constructor for creating a field attached to a mesh.
      * \param mesh PanNDE::Mesh* Raw pointer to the mesh (ownership transferred to field)
      * \param tp PanNDE::Field::FieldType The storage type of this field
      * \throw std::runtime_error If field type is not supported
      */
      HostField(PanNDE::Mesh* mesh,PanNDE::Field::FieldType tp){
        themesh.reset(mesh);
        setUp(tp);
      };
      
    private:
      /*!
      * Initializes internal data structures based on field type.
      * \param tp PanNDE::Field::FieldType The field type to set up
      * \throw std::runtime_error If field type is not implemented
      */
      void setUp(PanNDE::Field::FieldType tp){
        mytype=tp;
        field_data.resize(0);
        switch (mytype){
          case PanNDE::Field::NODE:
            field_data.resize(themesh->nodeCount());
            break;
          case PanNDE::Field::CELL:
            field_data.resize(themesh->cellCount());
            break;
          case PanNDE::Field::CONSTANT:
            field_data.resize(1);
            break;
          default:
            throw std::runtime_error("Field Type Not Implemented");
            break;
        };
      };
      
      //! The mesh to which this field is attached
      std::shared_ptr<PanNDE::Mesh> themesh;
      
      //! The storage type of this field (NODE, CELL, or CONSTANT)
      PanNDE::Field::FieldType mytype;
      
      //! The array storing field values
      std::vector<double> field_data;
  };

  /*! \class HostFieldBundle HostField.hpp "modules/HostData/include/HostField.hpp"
  *
  * Implements a container for a collection of named fields sharing a common mesh.
  * This class organizes multiple related fields and provides access by name
  * or index, maintaining their association with a parent mesh.
  *
  */
  class HostFieldBundle : public PanNDE::FieldBundle {
    public:
      /*!
      * Creates a shared pointer to a new empty HostFieldBundle instance.
      * \return std::shared_ptr<HostData::HostFieldBundle> Shared pointer to the new bundle
      */
      static
      std::shared_ptr<HostData::HostFieldBundle> makeShared(){
        auto bundle=std::make_shared<HostData::HostFieldBundle>(HostData::HostFieldBundle());
        return std::move(bundle);
      };
      
      /*!
      * Gets the parent mesh shared by all fields in this bundle.
      * \return std::shared_ptr<PanNDE::Mesh>& Reference to the mesh pointer
      */
      std::shared_ptr<PanNDE::Mesh>& mesh()override{return themesh;};
      
      /*!
      * Gets a field by its name.
      * \param keyname std::string The name of the field to retrieve
      * \return std::shared_ptr<PanNDE::Field> The requested field
      * \throw std::out_of_range If no field exists with the given name
      */
      std::shared_ptr<PanNDE::Field> field(std::string keyname)override{return fields.at(keyname);};
      
      /*!
      * Gets the name of a field by its index.
      * \param idx int The index of the field name to retrieve
      * \return std::string The name of the field at the specified index
      * \throw std::out_of_range If the index is out of bounds
      */
      std::string fieldName(int idx)override{return keynames.at(idx);};
      
      /*!
      * Gets the number of fields in this bundle.
      * \return int The number of fields
      */
      int fieldCount()override{return keynames.size();};
      
      /*!
      * Adds an existing field to the bundle with the given name.
      * \param keyname std::string Name to associate with the field
      * \param field std::shared_ptr<PanNDE::Field> The field to add
      * \throw std::runtime_error If a field with the given name already exists
      */
      void emplaceField(std::string keyname,std::shared_ptr<PanNDE::Field> field)override{
        auto result=fields.emplace(keyname,field);
        addKey(result.second,keyname);
      };
      
      /*!
      * Creates and adds a new field of the specified type to the bundle.
      * \param keyname std::string Name to associate with the new field
      * \param type PanNDE::Field::FieldType Type of field to create
      * \throw std::runtime_error If the bundle has no mesh or if a field with the given name already exists
      */
      void emplaceField(std::string keyname,PanNDE::Field::FieldType type)override{
        if(nullptr==themesh){throw std::runtime_error("Cannot construct field: No registered mesh with HostData::FieldBundle");};
        auto result=fields.emplace(keyname,
                            std::make_shared<HostData::HostField>(HostData::HostField(themesh,type)));
        addKey(result.second,keyname);
      };
      
      /*!
      * Default constructor creating an empty field bundle with no mesh.
      */
      HostFieldBundle(){};
      
      /*!
      * Constructor creating a field bundle associated with the specified mesh.
      * \param bundle_mesh std::shared_ptr<PanNDE::Mesh> The mesh for all fields in this bundle
      */
      HostFieldBundle(std::shared_ptr<PanNDE::Mesh> bundle_mesh){this->themesh=bundle_mesh;};
      
    private:
      /*!
      * Helper method to register a field name after insertion attempt.
      * \param inserted bool Whether the insertion was successful
      * \param keyname std::string The name that was used for insertion
      * \throw std::runtime_error If a field with the given name already exists
      */
      void addKey(bool inserted,std::string keyname){
        if(inserted){
          keynames.push_back(keyname);
        }else{throw std::runtime_error(keyname + " already in use.");};
      };
      
      //! The mesh shared by all fields in this bundle
      std::shared_ptr<PanNDE::Mesh> themesh=nullptr;
      
      //! Maps field names to their corresponding field objects
      std::map<std::string,std::shared_ptr<PanNDE::Field>> fields;
      
      //! Stores field names in insertion order
      std::vector<std::string> keynames;
  };

  /*! \class HostFieldFactory HostField.hpp "modules/HostData/include/HostField.hpp"
  * 
  * Factory class for creating host-based field objects and field bundles.
  * Implements the PanNDE::FieldFactory interface to provide methods for creating
  * various field-related objects.
  *
  */
  class HostFieldFactory : public PanNDE::FieldFactory {
    public:
      /*!
      * Creates a shared pointer to a new HostFieldFactory instance.
      * \return std::shared_ptr<HostData::HostFieldFactory> Shared pointer to the new factory
      */
      static
      std::shared_ptr<HostData::HostFieldFactory> makeShared(){
        auto maker=std::make_shared<HostData::HostFieldFactory>(HostData::HostFieldFactory());
        return std::move(maker);
      };
      
      /*!
      * Creates a new managed field attached to a mesh.
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to attach the field to
      * \param tp PanNDE::Field::FieldType The type of field to create
      * \return std::shared_ptr<PanNDE::Field> Shared pointer to the new field
      */
      std::shared_ptr<PanNDE::Field> makeManagedField(std::shared_ptr<PanNDE::Mesh> mesh,
                                                      PanNDE::Field::FieldType tp)override{
        auto field=std::make_shared<HostData::HostField>(HostData::HostField(mesh,tp));
        return std::move(field);
      };
      
      /*!
      * Creates a new managed field attached to a mesh.
      * \param mesh PanNDE::Mesh* Raw pointer to the mesh (ownership transferred to field)
      * \param tp PanNDE::Field::FieldType The type of field to create
      * \return std::shared_ptr<PanNDE::Field> Shared pointer to the new field
      */
      std::shared_ptr<PanNDE::Field> makeManagedField(PanNDE::Mesh* mesh,
                                                      PanNDE::Field::FieldType tp)override{
        auto field=std::make_shared<HostData::HostField>(HostData::HostField(mesh,tp));
        return std::move(field);
      };
      
      /*!
      * Creates a new raw field object (caller takes ownership).
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to attach the field to
      * \param tp PanNDE::Field::FieldType The type of field to create
      * \return PanNDE::Field* Pointer to the new field
      * \warning Use with caution - caller is responsible for memory management
      */
      PanNDE::Field* newField(std::shared_ptr<PanNDE::Mesh> mesh,
                              PanNDE::Field::FieldType tp)override{
        return new HostData::HostField(mesh,tp);
      };
      
      /*!
      * Creates a new raw field object (caller takes ownership).
      * \param mesh PanNDE::Mesh* Raw pointer to the mesh
      * \param tp PanNDE::Field::FieldType The type of field to create
      * \return PanNDE::Field* Pointer to the new field
      * \warning Use with caution - caller is responsible for memory management
      */
      PanNDE::Field* newField(PanNDE::Mesh* mesh,
                              PanNDE::Field::FieldType tp)override{
        return new HostData::HostField(mesh,tp);
      };

      /*!
      * Deletes a field created by newField().
      * \param field PanNDE::Field* Pointer to the field to delete
      */
      void deleteField(PanNDE::Field* field)override{if(nullptr!=field){delete field;};};

      /*!
      * Creates an empty array to store field objects.
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> Shared pointer to the new array
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> makeManagedFieldArray()override{
        return std::move(makeArray<std::shared_ptr<PanNDE::Field>>());
      };

      /*!
      * Creates an empty field bundle.
      * \return std::shared_ptr<PanNDE::FieldBundle> Shared pointer to the new field bundle
      */
      std::shared_ptr<PanNDE::FieldBundle> makeEmptyManagedFieldBundle()override{
        return std::move(std::make_shared<HostData::HostFieldBundle>(HostData::HostFieldBundle()));
      };
      
    private:
      /*!
      * Helper method to create an empty array.
      * \return std::shared_ptr<PanNDE::Array<T>> New empty array
      */
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> makeArray(){
        auto arraymfg=std::make_shared<HostData::HostArrayFactory<T>>(HostData::HostArrayFactory<T>());
        auto array=arraymfg->makeManagedArray();
        return std::move(array);
      };
  };
};