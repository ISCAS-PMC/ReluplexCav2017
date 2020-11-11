/*********************                                                        */
/*! \file absEle.h
 ** Top contributors (to current version):
 **   Jiangchao Liu
 ** This file is part of the AI & SMTReluplex project.
 ** Copyright (c) 2018-2100 by the authors listed in the file AUTHORS
 ** (in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **/

#ifndef __ABSELE_H
#define __ABSELE_H


#include <string>
#include <vector>
#include <iostream>
#include "apronxx/apronxx.hh"
#include <map>
#include "apronxx/apxx_box.hh"
#include "apronxx/apxx_oct.hh"
#include "apronxx/apxx_polka.hh"
#include "apronxx/apxx_t1p.hh"
#include "apronxx/apxx_texpr1.hh"
#include "apronxx/apxx_tcons1.hh"
#include "abst.h"

//#define ABST Abst_symV

static int unsplit = 0;
static int split = 0;
using namespace apron;

enum LogicOP {ABS_SUPEQ,ABS_SUP,ABS_EQ, ABS_DISEQ,ABS_LOWEQ,ABS_LOW};

class AbsEle {

 public:
  var genVar(int layer, int dim){
    return var(std::string("layer"
			   + std::to_string(layer)
			   + "dim" + std::to_string(dim)));
  }
  
  texpr1 vec2sub(int layer, int ndim, double *coffs){
    assert(ndim>0);
    environment env = _abst->get_environment();
    texpr1::builder firstCo(env,coffs[0]);
    texpr1::builder firstVar(env,genVar(layer,0));
    texpr1 expr(firstCo * firstVar);
    for (int j = 1; j < ndim; j++){
      expr = texpr1(texpr1::builder(expr) + 
		    texpr1::builder(env,coffs[j])
		    *texpr1::builder(env, genVar(layer,j)));
    }
    return texpr1(expr);
  }


  texpr1 vec2texpr(int layer, int ndim, double *coffs,  double consts){
    environment env = _abst->get_environment();
    texpr1 expr = vec2sub(layer, ndim, coffs);
    texpr1::builder bt = texpr1::builder(expr);
    return texpr1(bt+texpr1::builder(env,consts));
  }

   tcons1 vec2tcons(int layer, int ndim, double *coffs, 
		    LogicOP log, double consts){ 
     environment env = _abst->get_environment();
     texpr1 expr = vec2sub(layer, ndim, coffs);
     texpr1::builder bd(expr);
     texpr1::builder cst(env,consts);
     switch(log){
     case ABS_SUPEQ:
       return tcons1(expr >= cst);
     case ABS_LOWEQ:
       return tcons1(expr <= cst);
     case ABS_SUP:
       return tcons1(expr > cst);
     case ABS_EQ:
       return tcons1(expr == cst);
     case ABS_DISEQ:
       return tcons1(expr != cst);
     case ABS_LOW:
       return tcons1(expr < cst);
     default:
       throw "unexpected logic operator";
     }
  }
  
   AbsEle(ABST* aA){
     _abst = aA;
  }
  AbsEle(int ndim, std::vector<std::string> name, Domain dom, int num, double** coffs, LogicOP *lop, double **consts){
     _abst = new ABST(dom);
     for(int i = 0; i < ndim; i++) 
      _abst->add_var(genVar(0,i));
     for(int i = 0; i < num; i++)
      _abst->meet_tcons(vec2tcons(0, ndim,coffs[i],lop[i],*(consts[i])));
  }

  ~AbsEle(){
    delete _abst;
  }
  
  void matrix_assign(unsigned layer, int idim, int odim, double **weights, double **bias){
    for(unsigned i = 0; i < odim; i++){
      std::cout<<"layer"<<layer<<"dim"<<i<<std::endl;
      var dst = genVar(layer,i);
      _abst->add_var(dst);
      texpr1 expr = vec2texpr(layer-1, idim, weights[i], *bias[i]);
      _abst->assign(dst, expr);
      //   environment env = _abst->get_environment();
      //_abst->meet_tcons(tcons1(texpr1::builder(expr)==texpr1::builder(env, dst)));
    }    
    for(unsigned i = 0; i < idim; i++)
      _abst->rem_var(genVar(layer - 1, i));
  }
  
  void linRELU(unsigned layer, unsigned dim){
    var cvar = genVar(layer, dim);
    interval bound = _abst->get_bound(cvar);
    if(double(bound.get_sup()) <= 0 ){
      this->assign_zero(layer,dim);
      return;
    }
    if(double(bound.get_inf())>=0) return;    
    var tvar = var(std::string("tmp"));
    var ttvar = var(std::string("blayer"
			   + std::to_string(layer)
			   + "dim" + std::to_string(dim)));
    _abst->add_var(tvar);
    environment env = _abst->get_environment();
    texpr1::builder lh(env,double(bound.get_sup())*double(bound.get_inf()));
    texpr1::builder hMinusl(env,double(bound.get_sup())-double(bound.get_inf()));
    texpr1::builder h(env,double(bound.get_sup()));
    tcons1 consUp= h*texpr1::builder(env,cvar) - hMinusl*texpr1::builder(env,tvar)>=lh;
    tcons1 consDown= texpr1::builder(env,cvar) - texpr1::builder(env,tvar)<=texpr1::builder(env,0);
    tcons1 consPositive = texpr1::builder(env,tvar)>=texpr1::builder(env,0);
    _abst->meet_tcons(consUp);
    _abst->meet_tcons(consDown);
    _abst->meet_tcons(consPositive);
    _abst->rename(cvar,ttvar);
    _abst->rename(tvar,cvar);
  }

  
  

  void assert_positive(int layer, int dim){
    environment env = _abst->get_environment();
    texpr1::builder tVar(env,genVar(layer, dim));
    tcons1 tcons = tcons1(tVar >= texpr1::builder(env,0));
    _abst->meet_tcons(tcons);
    if (_abst->is_bottom()) unsplit+=2;
  }
  void assert_negative_strict(int layer, int dim){
    environment env = _abst->get_environment();
    texpr1::builder tVar(env, genVar(layer, dim));
    tcons1 tcons = tcons1(tVar < texpr1::builder(env, 0));
    _abst->meet_tcons(tcons);
    if (_abst->is_bottom()) unsplit+=2;
  }
  void assign_zero(int layer, int dim){
    environment env = _abst->get_environment();
    var dst = genVar(layer, dim);
    _abst->assign(dst, texpr1(texpr1::builder(env, 0)));
  }

  interval get_bound(unsigned layer, unsigned dim){
    var dst = genVar(layer, dim);
    return _abst->get_bound(dst);
  }
  
  std::vector<interval> get_bounds(unsigned layer, unsigned ndim){
    std::vector<interval> re;
    for(int i = 0; i < ndim; i++){
      interval tv = get_bound(layer, i);
      re.push_back(tv);
    }
    return re;
  }



  void join(AbsEle *abs){
    _abst->join(abs->_abst);
  }

  AbsEle* copy(){
    ABST *tabst = new ABST(_abst);
    return (new AbsEle(tabst)); 
  }

  friend std::ostream& operator<< (std::ostream& os, const AbsEle& a){
    std::cout<<*(a._abst)<<std::endl;
    return os;
  }

 private:
  ABST* _abst;
};

#endif
