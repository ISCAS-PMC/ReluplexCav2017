/*********************                                                        */
/*! \file absEngine.h
 ** Top contributors (to current version):
 **   Jiangchao Liu
 ** This file is part of the AI & SMTReluplex project.
 ** Copyright (c) 2018-2100 by the authors listed in the file AUTHORS
 ** (in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **/

#ifndef __ABSENGINE_H
#define __ABSENGINE_H


#include <string.h>
#include "dnn.h"
#include "absEle.h"
#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>

class AbsEngine{ 
 public:
 AbsEngine(Dnn *aDnn):_dnn(NULL) {
    _dnn = aDnn;
  }
  void iterateWholeNet(AbsEle * aAE){
    if (!_dnn || !aAE) throw "Engine falls";
    std::vector<Layer*> ly = _dnn->allLayers();
    for(unsigned i = 0; i < ly.size() ; i ++){
      //for(unsigned i = 0; i < 4 ; i ++){
      ly[i]->trans(i+1,aAE);
    }
  }
  Dnn* getDnn(){
    return _dnn;
 }
 ~AbsEngine() {
    delete _dnn;
  }
  
 private:
  Dnn *_dnn;
};
#endif
