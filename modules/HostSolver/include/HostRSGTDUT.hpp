/*! \headerfile HostRSGTDUT.hpp "modules/HostSolver/include/HostRSGTDUT.hpp"
* "HostRSGTDUT.hpp" contains the implementation of an elastic wave propagation model
* for ultrasonic testing using a Rotated Staggered Grid Time Domain (RSGTD) approach.
* This class provides an explicit, matrix-free finite difference solver for simulating
* ultrasonic wave propagation in arbitrary materials.
*
* Note: This implementation has been superseded by HostUTModel, which provides
* enhanced functionality including moment tensor excitation sources.
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
#include <array>
#include <set>

#include "StandardNames.hpp"

#include "Field.hpp"
#include "Mesh.hpp"
#include "MultiVariate.hpp"
#include "Communicator.hpp"
#include "Model.hpp"

#include "internal/RightHex.hpp"

namespace HostSolver {
  /*! \class HostRSGTDUT HostRSGTDUT.hpp "modules/HostSolver/include/HostRSGTDUT.hpp"
  * 
  * Implements an elastic wave propagation solver for ultrasonic testing using a 
  * rotated staggered-grid approach in the time domain.
  * 
  * This class provides a finite difference time domain solver for simulating 
  * elastic wave propagation in anisotropic materials. It uses a staggered-grid 
  * scheme where velocities are defined at nodes and stresses at cell centers, with
  * explicit time stepping. The solver supports parallel execution through domain
  * decomposition and halo exchanges.
  * 
  * \deprecated This class has been superseded by HostUTModel, which provides
  * additional capabilities including moment tensor excitation sources.
  *
  */
  class HostRSGTDUT : public PanNDE::Model {
    public:
      /*!
      * Creates a shared pointer to a new HostRSGTDUT instance.
      *
      * \param dt double Time step size
      * \param domain std::shared_ptr<PanNDE::FieldBundle> Bundle containing material properties
      * \param initial_condition std::shared_ptr<PanNDE::FieldBundle> Bundle with initial field values
      * \param transducers std::shared_ptr<PanNDE::MultiVariate>[][3] Array of transducer excitation functions
      * \param Nxducers int Number of transducers (default: 0)
      * \param comm std::shared_ptr<PanNDE::Communicator> Optional communicator for parallel execution
      * \return std::shared_ptr<HostSolver::HostRSGTDUT> Shared pointer to the new model
      */
      static std::shared_ptr<HostSolver::HostRSGTDUT> makeShared(double dt,
                  std::shared_ptr<PanNDE::FieldBundle> domain,
                  std::shared_ptr<PanNDE::FieldBundle> initial_condition,
                  std::shared_ptr<PanNDE::MultiVariate> transducers[][3],int Nxducers=0,
                  std::shared_ptr<PanNDE::Communicator> comm=nullptr){
      return std::move(std::make_shared<HostSolver::HostRSGTDUT>(
                HostSolver::HostRSGTDUT(dt,domain,initial_condition,
                                        transducers,Nxducers,comm)));
      };
      
      /*!
      * Constructs a model with a single transducer.
      *
      * \param dt double Time step size
      * \param domain std::shared_ptr<PanNDE::FieldBundle> Bundle containing material properties
      * \param initial_condition std::shared_ptr<PanNDE::FieldBundle> Bundle with initial field values
      * \param transducer std::shared_ptr<PanNDE::MultiVariate>[3] Array of force components for one transducer
      * \param comm std::shared_ptr<PanNDE::Communicator> Optional communicator for parallel execution
      */
      HostRSGTDUT(double dt,
                  std::shared_ptr<PanNDE::FieldBundle> domain,
                  std::shared_ptr<PanNDE::FieldBundle> initial_condition,
                  std::shared_ptr<PanNDE::MultiVariate> transducer[3],
                  std::shared_ptr<PanNDE::Communicator> comm=nullptr){
        setup(dt,domain,initial_condition,comm);
        resizeXDvector(1);
        addTransducer(0,transducer);
      };
      
      /*!
      * Constructs a model with multiple transducers.
      *
      * \param dt double Time step size
      * \param domain std::shared_ptr<PanNDE::FieldBundle> Bundle containing material properties
      * \param initial_condition std::shared_ptr<PanNDE::FieldBundle> Bundle with initial field values
      * \param transducers std::shared_ptr<PanNDE::MultiVariate>[][3] Array of transducer force components
      * \param Nxducers int Number of transducers (default: 0)
      * \param comm std::shared_ptr<PanNDE::Communicator> Optional communicator for parallel execution
      */
      HostRSGTDUT(double dt,
                  std::shared_ptr<PanNDE::FieldBundle> domain,
                  std::shared_ptr<PanNDE::FieldBundle> initial_condition,
                  std::shared_ptr<PanNDE::MultiVariate> transducers[][3],int Nxducers=0,
                  std::shared_ptr<PanNDE::Communicator> comm=nullptr){
        setup(dt,domain,initial_condition,comm);
        setupTransducers(transducers,Nxducers);
      };

      /*!
      * Gets the current solution fields (velocities and stresses).
      * \return std::shared_ptr<PanNDE::FieldBundle> Bundle containing current solution fields
      */
      std::shared_ptr<PanNDE::FieldBundle> getStates()override{return solution;};
      //std::shared_ptr<PanNDE::FieldBundle> getOutput()override{return output;};

      /*!
      * Advances the simulation by the specified number of time steps.
      * Each step updates velocities and then stresses, applying transducer forces.
      *
      * \param nsteps int Number of time steps to execute (default: 1)
      * \return double The new simulation time after stepping
      */
      double solve(int nsteps=1)override{
        for(int k=0;k<nsteps;k++){
          updateVs();
          updateTs();
          time+=dt;
        };
        return time;
      };

    private:
      /*!
      * Sets up the model with domain, initial conditions, and communication.
      *
      * \param dt double Time step size
      * \param domain std::shared_ptr<PanNDE::FieldBundle> Material properties bundle
      * \param initial_condition std::shared_ptr<PanNDE::FieldBundle> Initial conditions bundle
      * \param comm std::shared_ptr<PanNDE::Communicator> Optional parallel communicator
      */
      void setup(double dt,
                 std::shared_ptr<PanNDE::FieldBundle> domain,
                 std::shared_ptr<PanNDE::FieldBundle> initial_condition,
                 std::shared_ptr<PanNDE::Communicator> comm=nullptr){
        this->dt=dt;
        this->domain=domain;
        mesh=domain->mesh();

        setupProperties();
        setupSolution(initial_condition);
        setupComms(comm);
        ddqi=std::make_shared<HostSolver::RightHex>(HostSolver::RightHex(mesh));
      };
      
      /*!
      * Sets up material property fields from the domain bundle.
      * Extracts density and elasticity tensor components.
      */
      inline void setupProperties(){
        density=domain->field(property_names.density);
        for(int ki=0;ki<6;ki++){
          for(int kj=0;kj<6;kj++){
            CIJ[ki][kj]=domain->field(property_names.CIJ[ki][kj]);
          };
        };
      };
      
      /*!
      * Resizes the transducer vector to hold the specified number of transducers.
      * \param N int Number of transducers to accommodate
      */
      inline void resizeXDvector(int N=1){
        xducers.resize(0);xducers.resize(N);
      };
      
      /*!
      * Sets up the transducer array from the provided transducer data.
      *
      * \param transducers std::shared_ptr<PanNDE::MultiVariate>[][3] Array of transducer functions
      * \param Nxducers int Number of transducers to set up
      */
      inline void setupTransducers(std::shared_ptr<PanNDE::MultiVariate> transducers[][3],
                                   int Nxducers=0){
        resizeXDvector(Nxducers);
        for(int k=0;k<Nxducers;k++){
          addTransducer(k,transducers[k]);
        };
      };
      
      /*!
      * Sets up solution fields from the initial condition bundle.
      * Extracts velocity and stress field components.
      *
      * \param initial_condition std::shared_ptr<PanNDE::FieldBundle> Initial field values
      */
      inline void setupSolution(std::shared_ptr<PanNDE::FieldBundle> initial_condition){
        solution=initial_condition;
        for(int kv=0;kv<3;kv++){V[kv]=solution->field(state_names.V[kv]);};
        for(int ks=0;ks<6;ks++){S[ks]=solution->field(state_names.SI[ks]);};
      };
      
      /*!
      * Sets up communication for parallel execution.
      * Configures data exchange patterns for field values.
      *
      * \param comm std::shared_ptr<PanNDE::Communicator> Communicator for parallel execution
      */
      inline void setupComms(std::shared_ptr<PanNDE::Communicator> comm=nullptr){
        this->comm=comm;
        this->comm->setupDataLinks(solution);
        if(nullptr!=this->comm){getSendNodes();};
      };
      
      /*!
      * Adds a transducer to the specified position in the transducer array.
      *
      * \param idx int Index to place the transducer
      * \param transducer std::shared_ptr<PanNDE::MultiVariate>[3] Force components in x,y,z directions
      */
      inline void addTransducer(int idx,std::shared_ptr<PanNDE::MultiVariate> transducer[3]){
        for(int k=0;k<3;k++){xducers.at(idx).at(k)=transducer[k];};
      };

      /*!
      * Updates all stress components at all cells.
      */
      void updateTs(){
        for(int k=0;k<mesh->cellCount();k++){
          updateTsAtCell(k);
        };
      };
      
      /*!
      * Updates stress components at a specific cell.
      * Computes velocity gradients and applies constitutive law.
      *
      * \param c_idx int Cell index to update
      */
      void updateTsAtCell(int c_idx){
        double delV[3][3];
        for(int k1=0;k1<3;k1++){
          delV[0][k1]=ddqi->DDx_cc(c_idx,V[k1]);
          delV[1][k1]=ddqi->DDy_cc(c_idx,V[k1]);
          delV[2][k1]=ddqi->DDz_cc(c_idx,V[k1]);
        };
        for(int ks=0;ks<6;ks++){
          S[ks]->at(c_idx)+=dt*(CIJ[ks][0]->atCell(c_idx)*delV[0][0]+
                                CIJ[ks][1]->atCell(c_idx)*delV[1][1]+
                                CIJ[ks][2]->atCell(c_idx)*delV[2][2]+
                                CIJ[ks][3]->atCell(c_idx)*(delV[1][2]+delV[2][1])+
                                CIJ[ks][4]->atCell(c_idx)*(delV[0][2]+delV[2][0])+
                                CIJ[ks][5]->atCell(c_idx)*(delV[0][1]+delV[1][0]));
        };
      };

      /*!
      * Updates velocities, handling parallel communication if needed.
      */
      void updateVs(){
        if(nullptr==comm){
          updateAllVs();
        }else{
          //updateVsBySet(halonodes);
          updateAllVs();
          startExchange();
          //updateVsBySet(bulknodes);
          waitForExchange();
        };
      };
      
      /*!
      * Updates velocities at all nodes in the mesh.
      */
      void updateAllVs(){
        for(int k=0;k<mesh->nodeCount();k++){
          updateVsAtNode(k);
        };
      };
      
      /*!
      * Updates velocities for nodes in a specific set.
      * Used for handling bulk and halo nodes in parallel execution.
      *
      * \param nodeset std::set<int> Set of node indices to update
      */
      void updateVsBySet(std::set<int> nodeset){
        for(auto it=nodeset.begin();it!=nodeset.end();it++){
          //printf("[%i] node %i\n",comm->getProcessId(),*it);
          updateVsAtNode(*it);
        };
      };
      
      /*!
      * Updates velocity components at a specific node.
      * Computes stress divergence and applies it to velocities.
      *
      * \param n_idx int Node index to update
      */
      void updateVsAtNode(int n_idx){
        if(mesh->partitionId()==mesh->nodeHomePartition(n_idx)){
          double xd_args[4];
          mesh->nodeCoordinate(n_idx,xd_args);
          xd_args[3]=time;
          V[0]->at(n_idx)+=(dt/(density->atNode(n_idx)))*
                                                (ddqi->DDx_nc(n_idx,S[0])+
                                                 ddqi->DDy_nc(n_idx,S[5])+
                                                 ddqi->DDz_nc(n_idx,S[4])+
                                                 evalXds(0,n_idx,time));
                                                 //evalXds(0,xd_args));
          V[1]->at(n_idx)+=(dt/(density->atNode(n_idx)))*
                                                (ddqi->DDx_nc(n_idx,S[5])+
                                                 ddqi->DDy_nc(n_idx,S[1])+
                                                 ddqi->DDz_nc(n_idx,S[3])+
                                                 evalXds(1,n_idx,time));
                                                 //evalXds(1,xd_args));
          V[2]->at(n_idx)+=(dt/(density->atNode(n_idx)))*
                                                (ddqi->DDx_nc(n_idx,S[4])+
                                                 ddqi->DDy_nc(n_idx,S[3])+
                                                 ddqi->DDz_nc(n_idx,S[2])+
                                                 evalXds(2,n_idx,time));
                                                 //evalXds(2,xd_args));
        };
      };
      
      /*!
      * Evaluates transducer excitation at a point.
      * 
      * \param dir int Direction index (0-2)
      * \param xd_args double* Position and time (4D array)
      * \return double Total transducer force in the specified direction
      */
      inline double evalXds(int dir,double* xd_args){
        double result=0.;
        for(auto it=xducers.begin();it!=xducers.end();it++){
          result+=((nullptr==it->at(dir))?0.:(it->at(dir)->evalAt(xd_args)));
        };
        return result;
      };
      
      /*!
      * Evaluates transducer excitation at a node and time using cell centers.
      * 
      * \param dir int Direction index (0-2)
      * \param n_idx int Node index
      * \param time double Current simulation time
      * \return double Total transducer force in the specified direction
      */
      inline double evalXds(int dir,int n_idx,double time){
        double result=0.;
        int32_t cbox[8];
        mesh->connectedCells(n_idx,cbox);
        double xd_args[4];
        for(int kb=0;kb<8;kb++){
          if(-1!=cbox[kb]){
            cellcenter(cbox[kb],xd_args);
            xd_args[3]=time;
            for(auto it=xducers.begin();it!=xducers.end();it++){
              result+=((nullptr==it->at(dir))?0.:0.125*(it->at(dir)->evalAt(xd_args)));
            };
          };
        };
        return result;
      };
      /*inline bool checkXds(int dir,int n_idx){
        bool check=false;
        double coord[3];
        mesh->nodeCoordinate(n_idx,coord);
        for(auto it=xducers.begin();it!=xducers.end();it++){
          check=((nullptr==it->at(dir))?0.:0.125*(it->at(dir)->checkAt(xd_args)));
        };
        return check;
      };*/
      
      /*!
      * Computes the center coordinates of a cell.
      * 
      * \param c_idx int Cell index
      * \param cc double[3] Output array for cell center coordinates
      */
      inline void cellcenter(int c_idx,double cc[3]){
        int32_t box[8];
        double pt[3];
        cc[0]=0.;cc[1]=0.;cc[2]=0.;
        mesh->cell(c_idx,box);
        for(int kb=0;kb<8;kb++){
          mesh->nodeCoordinate(box[kb],pt);
          for(int kd=0;kd<3;kd++){cc[kd]+=0.125*pt[kd];};
        };
      };

      /*!
      * Identifies nodes that need to be sent to other processors.
      * Separates nodes into bulk (interior) and halo (boundary) sets.
      */
      void getSendNodes(){
        auto mesh=domain->mesh();
        int32_t box[8];
        bulknodes.clear();halonodes.clear();
        for(int k=0;k<mesh->cellCount();k++){
          mesh->cell(k,box);
          for(int kb=0;kb<8;kb++){
            if(mesh->partitionId()!=mesh->nodeHomePartition(box[kb])){
              for(int kbb=0;kbb<8;kbb++){halonodes.emplace(box[kbb]);};
              continue;
            };
          };
        };
        for(int k=0;k<mesh->nodeCount();k++){
          if(halonodes.end()!=halonodes.find(k)){bulknodes.emplace(k);};
        };
      };
      
      /*!
      * Initiates asynchronous exchange of velocity field values between processors.
      */
      void startExchange(){
        for(int kv=0;kv<3;kv++){comm->startHaloExchange(state_names.V[kv]);};
      };
      
      /*!
      * Waits for completion of velocity field exchange between processors.
      */
      void waitForExchange(){
        for(int kv=0;kv<3;kv++){comm->waitUntilDone(state_names.V[kv]);};
      };

      //! Standard names for elastic material properties
      PanNDE::ElasticMaterialNames property_names;
      
      //! Standard names for time domain metadata
      PanNDE::TimeDomainMetadataNames timing_names;
      
      //! Standard names for transducer excitation metadata
      PanNDE::TransducerExcitationNames xducer_names;
      
      //! Standard names for elastic state variables
      PanNDE::ElasticStateNames state_names;

      //! Time step size and current simulation time
      double dt,time=0.;
      
      //! Bundle containing material properties
      std::shared_ptr<PanNDE::FieldBundle> domain=nullptr;
      
      //! Mesh on which the simulation is performed
      std::shared_ptr<PanNDE::Mesh> mesh=nullptr;
      
      //! Elasticity tensor fields (Cijkl in Voigt notation)
      std::shared_ptr<PanNDE::Field> CIJ[6][6];
      
      //! Material density field
      std::shared_ptr<PanNDE::Field> density=nullptr;

      //! Bundle containing solution fields (velocities and stresses)
      std::shared_ptr<PanNDE::FieldBundle> solution=nullptr;
      
      //! Output bundle (commented out, not used)
      //std::shared_ptr<PanNDE::FieldBundle> output=nullptr;
      
      //! Velocity field components (Vx, Vy, Vz)
      std::shared_ptr<PanNDE::Field> V[3]={nullptr,nullptr,nullptr};
      
      //! Stress field components in Voigt notation (σxx, σyy, σzz, σyz, σxz, σxy)
      std::shared_ptr<PanNDE::Field> S[6]={nullptr,nullptr,nullptr,
                                           nullptr,nullptr,nullptr};

      //! Collection of transducer excitations (force components in x, y, z for each)
      std::vector<std::array<std::shared_ptr<PanNDE::MultiVariate>,3>> xducers;
      
      //! Spatial derivative operator
      std::shared_ptr<HostSolver::RightHex> ddqi=nullptr;
      
      //! Communicator for parallel execution
      std::shared_ptr<PanNDE::Communicator> comm=nullptr;
      
      //! Set of nodes at partition boundaries (halo)
      std::set<int> halonodes;
      
      //! Set of interior nodes (bulk)
      std::set<int> bulknodes;
  };
};