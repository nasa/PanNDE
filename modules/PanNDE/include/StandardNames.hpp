//Equivalent to the SIDS definitions to CGNS
//Add to this as more physics are implemented
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

#include <string>

namespace PanNDE {
  struct ElasticMaterialNames {
    /*
      1=xx
      2=yy
      3=zz
      4=yz=zy
      5=xz=zx
      6=xy=yx
    */
    std::string CIJ[6][6]={
                            {"C11","C12","C13","C14","C15","C16"},
                            {"C12","C22","C23","C24","C25","C26"},
                            {"C13","C23","C33","C34","C35","C36"},
                            {"C14","C24","C34","C44","C45","C46"},
                            {"C15","C25","C35","C45","C55","C56"},
                            {"C16","C26","C36","C46","C56","C66"}
                          };
    std::string Cijkl[3][3][3][3]={
                                    {//x
                                      {//x 
                                         //x     y     z
                                        {"C11","C16","C15"},//x
                                        {"C16","C12","C14"},//y
                                        {"C15","C14","C13"} //z
                                      },
                                      {//y 
                                         //x     y     z
                                        {"C16","C66","C56"},//x
                                        {"C66","C26","C46"},//y
                                        {"C56","C46","C36"} //z
                                      },
                                      {//z 
                                         //x     y     z
                                        {"C15","C56","C55"},//x
                                        {"C56","C25","C45"},//y
                                        {"C55","C45","C35"} //z
                                      },
                                    },
                                    {//y
                                      {//x 
                                         //x     y     z
                                        {"C16","C66","C56"},//x
                                        {"C66","C26","C46"},//y
                                        {"C56","C46","C36"} //z
                                      },
                                      {//y
                                         //x     y     z
                                        {"C12","C26","C25"},//x
                                        {"C26","C22","C24"},//y
                                        {"C25","C24","C23"} //z
                                      },
                                      {//z
                                         //x     y     z
                                        {"C14","C46","C45"},//x
                                        {"C46","C24","C44"},//y
                                        {"C45","C44","C34"} //z
                                      }
                                    },
                                    {//z
                                      {//x
                                         //x     y     z
                                        {"C15","C56","C55"},//x
                                        {"C56","C25","C45"},//y
                                        {"C55","C45","C35"} //z
                                      },
                                      {//y
                                         //x     y     z
                                        {"C14","C46","C45"},//x
                                        {"C46","C24","C44"},//y
                                        {"C45","C44","C34"} //z
                                      },
                                      {//z
                                         //x     y     z
                                        {"C13","C36","C35"},//x
                                        {"C36","C23","C34"},//y
                                        {"C35","C34","C33"} //z
                                      }
                                    }
                                 };
    std::string density="density";
  };
  struct ElasticStateNames {
    std::string V[3]={"Vx","Vy","Vz"};
    std::string S[3][3]={
                          {"Sxx","Sxy","Sxz"},
                          {"Sxy","Syy","Syz"},
                          {"Sxz","Syz","Szz"}
                        };
    std::string SI[6]={"Sxx","Syy","Szz","Syz","Sxz","Sxy"};
  };
  struct TransducerExcitationNames {
    std::string Coordinates[3]={"XD_X","XD_Y","XD_Z"};
    std::string TimeCoordinate="XD_T";
    std::string ForceValue[3]={"XD_Fx","XD_Fy","XD_Fz"};
  };
  struct TimeDomainMetadataNames {
    std::string dt="dt";
    std::string write_times="write_times";
  };


  class TransducerWindowedSineParameterNames {
    public:
      inline std::string transducerCount(){return "NTransducers";};

      inline std::string centerX(int xd_idx){return prefix+std::to_string(xd_idx)+"/XCenter";};
      inline std::string centerY(int xd_idx){return prefix+std::to_string(xd_idx)+"/YCenter";};
      inline std::string centerZ(int xd_idx){return prefix+std::to_string(xd_idx)+"/ZCenter";};
      inline std::string frequency(int xd_idx){return prefix+std::to_string(xd_idx)+"/Frequency";};
      inline std::string cycleCount(int xd_idx){return prefix+std::to_string(xd_idx)+"/NCycle";};
      inline std::string phase(int xd_idx){return prefix+std::to_string(xd_idx)+"/Phase";};
      inline std::string radius(int xd_idx){return prefix+std::to_string(xd_idx)+"/Radius";};
    private:
      std::string prefix="XD";
  };
};
