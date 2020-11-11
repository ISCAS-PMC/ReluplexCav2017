/*********************                                                        */
/*! \file dnn.h
 ** Top contributors (to current version):
 **   Jiangchao Liu
 ** This file is part of the AI & SMTReluplex project.
 ** Copyright (c) 2018-2100 by the authors listed in the file AUTHORS
 ** (in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **/

#ifndef __DNN_H
#define __DNN_H


#include <string.h>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "AcasNeuralNetwork.h"
#include "layer.h"

enum NetType {Acas, CNN};
enum LayerType {FN, CN, MAXP};


class Dnn{
 public:
 Dnn(NetType nt, const std::string &path): _dnet(NULL)
    {
      Log_info("BUild dnn");
      switch(nt)
	{
	case Acas:
	  Log_info("Acas network built in DNN");
	  _dnet = (void*)(new AcasNeuralNetwork(path));
	  break;
	case CNN:
	  throw "DNN initialization is called";
	  break;
	default:
	  throw "unknon Network in initialization";
	  break;
	}
      _tp = nt;
    }
  std::vector<Layer*> allLayers(){
    if (_tp == Acas){
      std::vector<Layer*> lys(((AcasNeuralNetwork *)_dnet)->getNumLayers());
      for(unsigned i = 0; i < lys.size(); i++){
	lys[i] = new FCLayer(i+1,((AcasNeuralNetwork *)_dnet)->getLayerSize(i),
		   ((AcasNeuralNetwork *)_dnet)->getLayerSize(i+1),
		   ((AcasNeuralNetwork *)_dnet)->getWeights(i),
		   ((AcasNeuralNetwork *)_dnet)->getBias(i));
      }
      return lys;
    }
    else
      throw "not implemented";
  }
  ~Dnn(){
    if (_tp == Acas) {
      AcasNeuralNetwork * adnet = (AcasNeuralNetwork *)_dnet;
      delete adnet;
    }
    else
      throw "not implemented";

  }
  
  int get_inputdim(){
    if (_tp == Acas) return ((AcasNeuralNetwork *)_dnet)->get_inputdim();
    else
      throw "not implemented";
  }
  double*  get_mins(){
    if (_tp == Acas){
     return  ((AcasNeuralNetwork *)_dnet)->getMins();
    }
    else
      throw "not implemented";
  }
  double*  get_maxes(){
    if (_tp == Acas){
      return ((AcasNeuralNetwork *)_dnet)->getMaxes();
    }
    else
      throw "not implemented";
  }


 private:
  void *_dnet;
  NetType _tp;
};
#endif
