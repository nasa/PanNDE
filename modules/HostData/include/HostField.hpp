/*! \headerfile HostField.hpp "modules/HostData/include/HostField.hpp"
* "HostField.hpp" contains the class implementation encapsulating 
* a data field attached to a mesh
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
  * Implements the expected characteristics of a field attached to a mesh. 
  *
  */
  class HostField : public PanNDE::Field {
    public:
      enum FieldType{
        NODE=0,
        CELL=1,
        CONSTANT=2//Note: This will be treated as cellcentered and restricted to nodes as req'd
      };
      /*!
      * get the type of field, i.e., how the data is stored (Cell/Nodal)
      */
      PanNDE::Field::FieldType type()override{return mytype;};
      /*!
      * get the mesh to which the field is attached
      */
      std::shared_ptr<PanNDE::Mesh> mesh()override{return themesh;};
      /*!
      * get the pointer to underlying data array. This should be used with caution, and only when
      * low level access is necessary (e.g., there is a need for a memcpy())
      */
      double* data()override{return field_data.data();};
      /*!
      * get the value by reference at the index. the index corresponds to the mesh index based on the 
      * whether the field is cell or node centered.
      *\param index location index of the value desired
      */
      double& at(int index)override{
        return ((PanNDE::Field::CONSTANT==mytype)?field_data.at(0):field_data.at(index));
      };

      /*!
      * get the value at cell index
      * \param cidx cell index desired
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
      * get the value at node index
      * \param nidx nodal index desired
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
      * get the number of values in the underlying data array
      */
      int32_t size()override{return field_data.size();};

      /*!
      * map data from other field into this field. This should copy if the other field is of the 
      * same type (cell/node), and this should map the data if the other is not the same type
      */
      void mapFrom(std::shared_ptr<PanNDE::Field> other)override{
        mapFrom(other.get());
      };
      /*!
      * map data from other field into this field. This should copy if the other field is of the 
      * same type (cell/node), and this should map the data if the other is not the same type.
      * This solution is not the preferred solution as it does not use a smart pointer.
      */
      void mapFrom(PanNDE::Field* other)override{
        themesh=other->mesh();
        if(PanNDE::Field::CONSTANT==mytype){setUp(other->type());};
        for(auto k=0;k<size();k++){
          at(k)=(PanNDE::Field::NODE==mytype)?other->atNode(k):other->atCell(k);
        };
      };

      /*!
      * constructor
      * \param mesh mesh on which the field is attached
      * \param tp field type
      */
      HostField(std::shared_ptr<PanNDE::Mesh> mesh,PanNDE::Field::FieldType tp){
        themesh=mesh;
        setUp(tp);
      };
      /*!
      * constructor
      * \param mesh mesh on which the field is attached
      * \param tp field type
      */
      HostField(PanNDE::Mesh* mesh,PanNDE::Field::FieldType tp){
        themesh.reset(mesh);
        setUp(tp);
      };
    private:
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
      std::shared_ptr<PanNDE::Mesh> themesh;
      PanNDE::Field::FieldType mytype;
      std::vector<double> field_data;
  };

  /*! \class HostFieldBundle HostField.hpp "modules/HostData/include/HostField.hpp"
  *
  * Implements a container for passing related fields with their parent mesh
  *
  */
  class HostFieldBundle : public PanNDE::FieldBundle {
    public:
      /*!
      * get the parent mesh
      */
      std::shared_ptr<PanNDE::Mesh>& mesh()override{return themesh;};
      /*!
      * get the field by name
      * \param keyname field name
      */
      std::shared_ptr<PanNDE::Field> field(std::string keyname)override{return fields.at(keyname);};
      /*!
      * get the field name by index
      * \param idx index of the field in underlying storage
      */
      std::string fieldName(int idx)override{return keynames.at(idx);};
      /*!
      * get the number of bundled fields
      */
      int fieldCount()override{return keynames.size();};
      /*!
      * emplace extant field
      * \param keyname field name
      * \param field field to be emplaced
      */
      void emplaceField(std::string keyname,std::shared_ptr<PanNDE::Field> field)override{
        auto result=fields.emplace(keyname,field);
        addKey(result.second,keyname);
      };
      /*!
      * emplace new field
      * \param keyname field name
      * \param type field type to be constructed
      */
      void emplaceField(std::string keyname,PanNDE::Field::FieldType type)override{
        if(nullptr==themesh){throw std::runtime_error("Cannot construct field: No registered mesh with HostData::FieldBundle");};
        auto result=fields.emplace(keyname,
                            std::make_shared<HostData::HostField>(HostData::HostField(themesh,type)));
        addKey(result.second,keyname);
      };
      /*!
      * default constructor
      */
      HostFieldBundle(){};
      /*!
      * constructor
      * \param mesh for all fields
      */
      HostFieldBundle(std::shared_ptr<PanNDE::Mesh> bundle_mesh){this->themesh=bundle_mesh;};
    private:
      void addKey(bool inserted,std::string keyname){
        if(inserted){
          keynames.push_back(keyname);
        }else{throw std::runtime_error(keyname + " already in use.");};
      };
      std::shared_ptr<PanNDE::Mesh> themesh=nullptr;
      std::map<std::string,std::shared_ptr<PanNDE::Field>> fields;
      std::vector<std::string> keynames;
  };

  /*! \class HostFieldFactory HostField.hpp "modules/HostData/include/HostField.hpp"
  * Implements a factory class to create the HostField class
  *
  */
  class HostFieldFactory : public PanNDE::FieldFactory {
    public:
      /*! 
      * create a shared field
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      std::shared_ptr<PanNDE::Field> makeManagedField(std::shared_ptr<PanNDE::Mesh> mesh,
                                                      PanNDE::Field::FieldType tp)override{
        auto field=std::make_shared<HostData::HostField>(HostData::HostField(mesh,tp));
        return std::move(field);
      };
      /*! 
      * create a shared field
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      std::shared_ptr<PanNDE::Field> makeManagedField(PanNDE::Mesh* mesh,
                                                      PanNDE::Field::FieldType tp)override{
        auto field=std::make_shared<HostData::HostField>(HostData::HostField(mesh,tp));
        return std::move(field);
      };
      /*! 
      * create a field
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      PanNDE::Field* newField(std::shared_ptr<PanNDE::Mesh> mesh,
                              PanNDE::Field::FieldType tp)override{
        return new HostData::HostField(mesh,tp);
      };
      /*! 
      * create a field
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      PanNDE::Field* newField(PanNDE::Mesh* mesh,
                              PanNDE::Field::FieldType tp)override{
        return new HostData::HostField(mesh,tp);
      };

      /*! 
      * delete a field created with newField(mesh,tp)
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      void deleteField(PanNDE::Field* field)override{if(nullptr!=field){delete field;};};

      /*! 
      * create an array of fields
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> makeManagedFieldArray()override{
        return std::move(makeArray<std::shared_ptr<PanNDE::Field>>());
      };

      /*! 
      * create a field bundle
      */
      std::shared_ptr<PanNDE::FieldBundle> makeEmptyManagedFieldBundle(){
        return std::move(std::make_shared<HostData::HostFieldBundle>(HostData::HostFieldBundle()));
      };
    private:
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> makeArray(){
        auto arraymfg=std::make_shared<HostData::HostArrayFactory<T>>(HostData::HostArrayFactory<T>());
        auto array=arraymfg->makeManagedArray();
        return std::move(array);
      };
  };
};