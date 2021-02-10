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
#include <memory>
#include <string>
#include <cstdio>
#include <stdexcept>
#include <fstream>
#include <cstring>

namespace ParamUtils {
  struct ParsedLine{
    std::string name="";
    std::string content="";
    std::vector<double> fdata;
    std::vector<int> idata;
    int start_idx=-1;
  }; 
  class Node {
    public:
      static std::shared_ptr<Node> MakeShared(std::shared_ptr<Node> parent,
                                              ParsedLine node_data){
        return std::make_shared<Node>(Node(parent,node_data));
      };
      Node(std::shared_ptr<Node> parent,ParsedLine node_data){
        this->parent=parent;
        name=node_data.name;
        content=node_data;
        children.resize(0);
      };

      void assignChild(std::shared_ptr<Node> child){children.push_back(child);};
      std::shared_ptr<Node> getParent(){return parent;};
      std::shared_ptr<Node> getChild(int idx){return children.at(idx);};

      std::shared_ptr<Node> FindChildByName(std::string child_name){
        std::shared_ptr<Node> child;
        int child_idx=0;
        auto Nchildren=getNumberOfChildren();
        size_t sidx=std::string::npos;
        while(std::string::npos==sidx && child_idx<Nchildren){
          child=getChild(child_idx);
          sidx=child->getName().find(child_name);
          child_idx++;
        };
        if(std::string::npos==sidx){throw std::runtime_error("no child named "+child_name+" found");};
        return child;
      };

      size_t getNumberOfChildren(){return children.size();};

      std::string getName(){return name;};
      ParsedLine getContents(){return content;};

      double GetScalarD(){
        return getContents().fdata.at(0);
      };
      int GetScalarI(){
        return getContents().idata.at(0);
      };

    private:
      std::string name;
      ParsedLine content;
      std::vector<std::shared_ptr<Node>> children;
      std::shared_ptr<Node> parent;
  };
  class YAMLReader {
    public:
      static std::shared_ptr<YAMLReader> MakeShared(){
        return std::make_shared<YAMLReader>(YAMLReader());
      };

      int read(std::string filename){
        std::ifstream infile(filename);
        std::string line;

        while(!infile.eof()){
          std::getline(infile,line);
          if(0!=strcmp("#",line.substr(0,1).c_str()) && 0!=line.size()){
            contents.push_back(line);
          };
        };
        return lineCount();
      };
      int lineCount(){return contents.size();};
      std::string getLine(int idx){return contents.at(idx);};
    private:
      std::vector<std::string> contents;
  };
  class Parser {
    public:
      static std::shared_ptr<Parser> MakeShared(){
        return std::make_shared<Parser>(Parser());
      };

      
      ParsedLine parse(std::string raw_line){
        ParsedLine line_data;
        line_data.start_idx=raw_line.find_first_not_of(" ");
        line_data.name=getEntryName(raw_line);
        line_data.content=getDataString(raw_line);
        line_data.fdata=getData<double>(raw_line);
        line_data.idata=getData<int>(raw_line);
        return std::move(line_data);
      };

    private:
      size_t findSeparator(std::string raw_line){
        auto idx1=raw_line.find(":");
        if(raw_line.npos==idx1){throw std::runtime_error("line must have \":\"");};
        if(raw_line.find(": ")!=idx1 && idx1+1<raw_line.size()){
          throw std::runtime_error("colon must be followed by space");
        };
        return idx1;
      };
      std::string trim(std::string line,int start,int end){
        return line.substr(start,end-start);
      };
      std::string trimWhitespace(std::string line){
        auto sidx=line.find_first_not_of(" ");
        auto eidx=line.find_last_not_of(" ");
        return trim(line,sidx,eidx+1);
      };
      std::string getEntryName(std::string line){
        auto sep_idx=findSeparator(line);
        auto start_idx=line.find_first_not_of(" ");
        return trim(line,start_idx,sep_idx);
      };
      std::string getDataString(std::string line){
        std::string datastr="";
        auto sep_idx=findSeparator(line);
        if(line.size()>(sep_idx+2)){datastr=trimWhitespace(line.substr(sep_idx+2));};
        return datastr;
      };
      bool isVectorData(std::string content){
        auto br1=content.find_first_of("[");
        auto br2=content.find_last_of("]");
        return (0==br1 && content.size()-1==br2);
      };
      template<typename T>
      std::vector<T> getData(std::string raw_line){
        std::string content=getDataString(raw_line);
        if(isVectorData(content)){return parseVector<T>(content);};
        return parseScalar<T>(content);
      };
      template<typename T>
      std::vector<T> parseScalar(std::string content){
        std::vector<T> data;
        T value;
        if(0==castToNumber<T>(&value,content)){data.push_back(value);};
        return data;
      };
      template<typename T>
      std::vector<T> parseVector(std::string content){
        auto sdata=trim(content,1,content.size()-1);
        size_t sep_idx=0;
        std::vector<T> data;
        T value;
        do{
          sep_idx=sdata.find_first_of(",");
          auto entry=trim(sdata,0,sep_idx);
          if(0==castToNumber<T>(&value,entry)){data.push_back(value);};
          sdata=trim(sdata,sep_idx+1,sdata.size());
        }while(sdata.npos!=sep_idx);
        return data;
      };
      template<typename T>
      int castToNumber(T* value,std::string entry){
        throw std::runtime_error("invalid type");
        return 1;
      };
  };
  template<>
  int Parser::castToNumber<double>(double* value,std::string entry){
    try{
      *value=std::stod(entry);
      return 0;
    }catch(...){return 1;};
  };
  template<>
  int Parser::castToNumber<int>(int* value,std::string entry){
    try{
      *value=std::stoi(entry);
      return 0;
    }catch(...){return 1;};
  };


  std::vector<std::shared_ptr<Node>> BuildTreeFromFile(std::string filename){
    std::vector<std::shared_ptr<Node>> tree;
    auto reader=YAMLReader::MakeShared();
    auto parser=Parser::MakeShared();
    reader->read(filename);
    auto Nlines=reader->lineCount();
    ParsedLine root_node;root_node.name="ROOT";
    tree.push_back(Node::MakeShared(nullptr,root_node));
    std::vector<int> heritage_indent_level;
    std::vector<std::shared_ptr<Node>> current_heritage;
    current_heritage.resize(1);
    current_heritage.at(0)=tree.back();
    heritage_indent_level.resize(1);
    heritage_indent_level.at(0)=-2;

    for(int k=0;k<Nlines;k++){
      auto node_string=reader->getLine(k);
      auto node_data=parser->parse(node_string);
      int indent_level=node_data.start_idx;
      while(indent_level<=heritage_indent_level.back()){
        current_heritage.erase(current_heritage.end()-1);
        heritage_indent_level.erase(heritage_indent_level.end()-1);
      };
      tree.push_back(Node::MakeShared(current_heritage.back(),node_data));
      current_heritage.back()->assignChild(tree.back());
      heritage_indent_level.push_back(indent_level);
      current_heritage.push_back(tree.back());
    };
    return std::move(tree);
  };

  void PrintTree(std::vector<std::shared_ptr<Node>> tree){
    for(int node_idx=0;node_idx<tree.size();node_idx++){
      std::vector<std::shared_ptr<Node>> heritage;
      heritage.resize(0);heritage.push_back(tree.at(node_idx)->getParent());
      while(nullptr!=heritage.back()){
        heritage.push_back(heritage.back()->getParent());
      };
      for(int k=1;k<heritage.size();k++){
        printf("\'%s\'\n  |\n",heritage.at(heritage.size()-1-k)->getName().c_str());
      };
      printf("\'%s\'\n",tree.at(node_idx)->getName().c_str());
      if(0<tree.at(node_idx)->getNumberOfChildren()){printf("|||||||\n");};
      std::string child_names;
      for(int kc=0;kc<(tree.at(node_idx)->getNumberOfChildren());kc++){
        child_names+="\'"+tree.at(node_idx)->getChild(kc)->getName()+"\'  ";
      };
      printf("%s\n\n",child_names.c_str());
    };
  };
};