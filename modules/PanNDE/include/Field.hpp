/*! \headerfile Field.hpp "modules/PanNDE/include/Field.hpp"
* "Field.hpp" contains the class definition encapsulating 
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

#include <cstdint>
#include <memory>

#include "Mesh.hpp"

namespace PanNDE{
  /*! \class Field Field.hpp "modules/PanNDE/include/Field.hpp"
  *
  * Defines the expected characteristics of a field attached to a mesh. 
  *
  */
  class Field{
    public:
      enum FieldType{
        NODE=0,
        CELL=1,
        CONSTANT=2//Note: This will be treated as cellcentered and restricted to nodes as req'd
      };
      /*!
      * get the type of field, i.e., how the data is stored (Cell/Nodal)
      */
      virtual FieldType type() =0;
      /*!
      * get the mesh to which the field is attached
      */
      virtual std::shared_ptr<PanNDE::Mesh> mesh() =0;
      /*!
      * get the pointer to underlying data array. This should be used with caution, and only when
      * low level access is necessary (e.g., there is a need for a memcpy())
      */
      virtual double* data() =0;
      /*!
      * get the value by reference at the index. the index corresponds to the mesh index based on the 
      * whether the field is cell or node centered.
      *\param index location index of the value desired
      */
      virtual double& at(int index) =0;

      /*!
      * get the value at cell index
      * \param cidx cell index desired
      */
      virtual double atCell(int cidx) =0;
      /*!
      * get the value at node index
      * \param nidx nodal index desired
      */
      virtual double atNode(int nidx) =0;

      /*!
      * get the number of values in the underlying data array
      */
      virtual int32_t size() =0;

      /*!
      * map data from other field into this field. This should copy if the other field is of the 
      * same type (cell/node), and this should map the data if the other is not the same type
      */
      virtual void mapFrom(std::shared_ptr<PanNDE::Field> other) =0;
      /*!
      * map data from other field into this field. This should copy if the other field is of the 
      * same type (cell/node), and this should map the data if the other is not the same type.
      * This solution is not the preferred solution as it does not use a smart pointer.
      */
      virtual void mapFrom(PanNDE::Field* other) =0;
  };

  /*! \class FieldBundle Field.hpp "modules/PanNDE/include/Field.hpp"
  *
  * Defines a container for passing related fields with their parent mesh
  *
  */
  class FieldBundle {
    public:
      /*!
      * get the parent mesh
      */
      virtual std::shared_ptr<PanNDE::Mesh>& mesh() =0;
      /*!
      * get the field by name
      * \param keyname field name
      */
      virtual std::shared_ptr<PanNDE::Field> field(std::string keyname) =0;
      /*!
      * get the field name by index
      * \param idx index of the field in underlying storage
      */
      virtual std::string fieldName(int idx) =0;
      /*!
      * get the number of bundled fields
      */
      virtual int fieldCount() =0;
      /*!
      * emplace extant field
      * \param keyname field name
      * \param field field to be emplaced
      */
      virtual void emplaceField(std::string keyname,std::shared_ptr<PanNDE::Field> field) =0;
      /*!
      * emplace new field
      * \param keyname field name
      * \param type field type to be constructed
      */
      virtual void emplaceField(std::string keyname,PanNDE::Field::FieldType type) =0;
  };

  /*! \class FieldFactory Field.hpp "modules/PanNDE/include/Field.hpp"
  *
  * Defines a factory class to create the Field class
  *
  */
  class FieldFactory{
    public:
      /*! 
      * create a shared field
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      virtual std::shared_ptr<PanNDE::Field> makeManagedField(std::shared_ptr<PanNDE::Mesh> mesh,
                                                              PanNDE::Field::FieldType tp) =0;
      /*! 
      * create a shared field
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      virtual std::shared_ptr<PanNDE::Field> makeManagedField(PanNDE::Mesh* mesh,
                                                              PanNDE::Field::FieldType tp) =0;
      /*! 
      * create a field
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      virtual PanNDE::Field* newField(std::shared_ptr<PanNDE::Mesh> mesh,
                                      PanNDE::Field::FieldType tp) =0;
      /*! 
      * create a field
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      virtual PanNDE::Field* newField(PanNDE::Mesh* mesh,
                                      PanNDE::Field::FieldType tp) =0;
      /*! 
      * delete a field created with newField(mesh,tp)
      * \param mesh mesh on which field is constructed
      * \param tp type of field to make
      */
      virtual void deleteField(PanNDE::Field* field) =0;
      
      /*! 
      * create an array of fields
      */
      virtual std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> makeManagedFieldArray() =0;

      /*! 
      * create a field bundle
      */
      virtual std::shared_ptr<PanNDE::FieldBundle> makeEmptyManagedFieldBundle() =0;
  };
};


