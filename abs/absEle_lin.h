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
#include "Log.h"
#include "apronxx/apxx_box.hh"
#include "apronxx/apxx_oct.hh"
#include "apronxx/apxx_polka.hh"
#include "apronxx/apxx_t1p.hh"
#include "apronxx/apxx_texpr1.hh"
#include "apronxx/apxx_tcons1.hh"

static int unsplit = 0;
static int split = 0;
using namespace apron;

enum Domain {BOX, OCT, POLKA, TOPO};
enum LogicOP {ABS_SUPEQ,ABS_SUP,ABS_EQ, ABS_DISEQ,ABS_LOWEQ,ABS_LOW};

std::map<LogicOP,ap_constyp_t> logAp{ {ABS_SUPEQ,AP_CONS_SUPEQ},
    {ABS_SUP,AP_CONS_SUP}, {ABS_EQ,AP_CONS_EQ},
    {ABS_DISEQ,AP_CONS_DISEQ}, {ABS_LOWEQ,AP_CONS_SUPEQ},
    {ABS_LOW,AP_CONS_SUP}
};

std::map<LogicOP,int> logSign{ {ABS_SUPEQ,1},
    {ABS_SUP,1}, {ABS_EQ,1},
    {ABS_DISEQ,1}, {ABS_LOWEQ,-1},
    {ABS_LOW,-1}
};

class AbsEle {

 public:

  linexpr1 vec2expr(int layer, int ndim, double *coffs,  double consts){
    environment env = _elt->get_environment();
    coeff *tcoeff = new coeff[ndim];		
    unsigned realDim;
    for (int j = 0; j < ndim; j++){
      var cur = var(std::string("layer" + 
			     std::to_string(layer)
			     +
			     "dim"+std::to_string(j)));
      realDim = env.get_dim(cur);
       tcoeff[realDim] = coeff(coffs[j]);
     }
     linexpr1  lexpr1(env,linexpr0(ndim,tcoeff,coeff(consts)));
     delete[] tcoeff;
     return lexpr1;
   }

  
   std::vector<const linexpr1*>  matrix2lexpr1(int layer, int ndim, int num,
				 double** coffs,
				 double **consts){ 
    std::vector<const linexpr1*> t_lexpr(num);
    for(int i = 0; i < num; i++){
      linexpr1 lexpr1 = vec2expr(layer, ndim,coffs[i],*(consts[i]));
      t_lexpr[i] = new linexpr1(lexpr1);
    }
    return t_lexpr;
  }
  
  
  
   lincons1 vec2lcon(int layer, int ndim, double *coffs, 
		    LogicOP log, double consts){ 
    environment env = _elt->get_environment();
    coeff *tcoeff = new coeff[ndim];
    int sgn = logSign[log];
    unsigned realDim;
    for (int j = 0; j < ndim; j++){
      var cur = var(std::string("layer"
			    + std::to_string(layer)
			    + "dim" + std::to_string(j)));
      realDim = env.get_dim(cur);
      tcoeff[realDim] = coeff(sgn*coffs[j]);
    }
    linexpr0 t_lexpr0(ndim,tcoeff,coeff(-sgn*(consts)));
    lincons0 t_lcons0(logAp[log],t_lexpr0);
    delete[] tcoeff;
    lincons1 lcon1(env, t_lcons0);
    return lcon1;
  }
  

   lincons1_array matrix2lcons(int layer, int ndim, int num,
			      double** coffs, LogicOP *lop,
			      double **consts){
    std::vector<lincons1> lcons1Array;
    for(int i = 0; i < num; i++)
      lcons1Array.push_back(vec2lcon(layer, ndim,coffs[i],lop[i],*(consts[i])));
    return lincons1_array(lcons1Array);
  }

   tcons1_array matrix2tcons(int layer, int ndim, int num,
			    double** coffs, LogicOP *lop,
			       double **consts){
     std::vector<tcons1> tcons1Array;
     for(int i = 0; i < num; i++)
       tcons1Array.push_back(tcons1(vec2lcon(layer, ndim,coffs[i],lop[i],*(consts[i]))));
     return tcons1_array(tcons1Array);
   }
  
  
  AbsEle( manager *man, std::vector<var> dims, abstract1 *elt){
    _dims = dims;
    _elt = elt;
    _man = man;
  }

  /* AbsEle(int ndim, std::vector<std::string> name, Domain dom, int num, double** coffs, LogicOP *lop, double **consts){ */
  /*   //set environment */
  /*   Log_info("set environment");  */
  /*   for(int i = 0; i < ndim; i++)  */
  /*     _dims.push_back(var(std::string("layer0")+ "dim" + std::to_string(i))); */
  /*   std::vector<var> t_intdims; */
  /*   environment env = environment(t_intdims,_dims); */
  /*   Log_info("set manager"); */
  /*   // set manager */
  /*   switch (dom){ */
  /*   case BOX: */
  /*     _man = new box_manager(); */
  /*     break; */
  /*   case OCT: */
  /*     _man = new oct_manager(); */
  /*     break; */
  /*   case POLKA: */
  /*     _man = new polka_manager(); */
  /*     break; */
  /*   case TOPO: */
  /*     _man = new t1p_manager(); */
  /*     break; */
  /*   default: */
  /*     throw "not implemented"; */
  /*   } */
  /*   // build abstract ele with  constraints */
  /*   _elt = new abstract1(*_man,env,top()); */
  /*   tcons1_array tcons1Array =matrix2tcons(0, ndim, num,coffs,lop,consts); */
  /*   std::cout<<*_elt<<std::endl; */
  /*   std::cout<<tcons1Array<<std::endl; */
  /*   _elt->meet(*_man,tcons1Array); */
  /*   std::cout<<"after initialization \n"<<*_elt<<std::endl; */
  /* } */


  AbsEle(int ndim, std::vector<std::string> name, Domain dom, int num, double** coffs, LogicOP *lop, double **consts){
    //set environment
    Log_info("set environment"); 
    for(int i = 0; i < ndim; i++) 
      _dims.push_back(var(std::string("layer0")+ "dim" + std::to_string(i)));
    std::vector<var> t_intdims;
    environment env = environment(t_intdims,_dims);
    Log_info("set manager");
    // set manager
    switch (dom){
    case BOX:
      _man = new box_manager();
      break;
    case OCT:
      _man = new oct_manager();
      break;
    case POLKA:
      _man = new polka_manager();
      break;
    case TOPO:
      _man = new t1p_manager();
      break;
    default:
      throw "not implemented";
    }
    // build abstract ele with  constraints
    _elt = new abstract1(*_man,env,top());
    tcons1_array tcons1Array =matrix2tcons(0, ndim, num,coffs,lop,consts);
    std::cout<<*_elt<<std::endl;
    std::cout<<tcons1Array<<std::endl;
    _elt->meet(*_man,tcons1Array);
    std::cout<<"after initialization \n"<<*_elt<<std::endl;
  }

  ~AbsEle(){
    delete _elt;
    delete _man;
  }
  
  void matrix_assign(unsigned layer, int idim, int odim, double **weights, double **bias){
    std::vector<var> dstdims;
    environment env = _elt->get_environment();
    //std::cout<<"dstdims"<<std::endl;
    for(int i = 0; i < odim; i++)
      dstdims.push_back(var(std::string("layer"
					+ std::to_string(layer)
					+ "dim" + std::to_string(i))));
    std::vector<var> tt;    
    env = env.add(tt,dstdims);
    _elt->change_environment(*_man,env);
    std::vector<const linexpr1*> lexpr = 
      matrix2lexpr1(layer-1,idim, odim,weights,bias);
    for(int i = 0; i< odim; i++){ 
      var dst = var(std::string("layer"
			    + std::to_string(layer)
			    + "dim" + std::to_string(i)));
      _elt->assign(*_man,dst,texpr1(*(lexpr[i])));
    }
    env = env.remove(_dims);
    _dims = dstdims;
    for(int i = 0; i< odim; i++) delete lexpr[i];
     _elt->change_environment(*_man,env);
  }
  
  void linRELU(unsigned layer, unsigned dim){
    std::cout<<"layer "<<layer<<" dim "<< dim<<std::endl;
    var cvar = var(std::string("layer"
			       + std::to_string(layer)
			       + "dim" + std::to_string(dim)));
    interval bound = _elt->bound(*_man,cvar);
    if(double(bound.get_sup()) <= 0 ){
      this->assign_zero(layer,dim);
      return;
    }
    if(double(bound.get_inf())>=0) return;
    environment env = _elt->get_environment();
    var tvar = var(std::string("tmp"));
    std::vector<var> ti;    
    std::vector<var> tf;    
    tf.push_back(tvar);
    env = env.add(ti,tf);
    _elt->change_environment(*_man,env);
    texpr1::builder lh(env,double(bound.get_sup())*double(bound.get_inf()));
    texpr1::builder hMinusl = texpr1::builder(env,double(bound.get_sup())-double(bound.get_inf()));
    texpr1::builder h = texpr1::builder(env,double(bound.get_sup()));
    tcons1 consUp= h*texpr1::builder(env,cvar) - hMinusl*texpr1::builder(env,tvar)>=lh;
    tcons1 consDown= texpr1::builder(env,cvar) - texpr1::builder(env,tvar)<=texpr1::builder(env,0);
    tcons1 consPositive = texpr1::builder(env,tvar)>=texpr1::builder(env,0);
    std::vector<tcons1> tconsVec = {consUp, consDown, consPositive};
    tcons1_array tarray = tcons1_array(tconsVec);
    _elt->meet(*_man,tarray);
    std::vector<var> tc;
    tc.push_back(cvar);
    env = env.remove(tc);
    _elt->change_environment(*_man,env);
    _elt->rename(*_man,tf,tc);
  }



  void assert_positive(int layer, int dim){
    //std::cout<<"ASSERT:+ begin "<<std::endl;
    double *coff = new double[_dims.size()];
    memset(coff, 0, _dims.size() * sizeof(double));
    coff[dim] = 1;
    LogicOP log = ABS_SUPEQ;
    std::vector<tcons1> tconsVec = {tcons1(vec2lcon(layer, _dims.size(), coff, log, 0))};
    tcons1_array tarray = tcons1_array(tconsVec);
    _elt->meet(*_man,tarray);
    if (_elt->is_bottom(*_man)) unsplit++;
    // else split++;
    // Li Jianlin

    //std::cout<<"after is bottom"<<std::endl;
    delete[] coff;
  }
  void assert_negative_strict(int layer, int dim){
    //std::cout<<"ASSERT:- begin "<<std::endl;
    double *coff = new double[_dims.size()];
    memset(coff,0,_dims.size() * sizeof(double));
    coff[dim] = 1;
    LogicOP log = ABS_LOW;
    std::vector<tcons1> tconsVec = {tcons1(vec2lcon(layer, _dims.size(), coff, log, 0))};
    tcons1_array tarray = tcons1_array(tconsVec);
    //std::cout<<"ASSERT:- "<<larray<<std::endl;
    _elt->meet(*_man,tarray);
    if (_elt->is_bottom(*_man)) unsplit++;
    // else split++;
    // Li Jianlin

    //std::cout<<"after meet"<<std::endl;
    delete[] coff;
    //std::cout<<"delete"<<std::endl;
  }
  void assign_zero(int layer, int dim){
    environment env = _elt->get_environment();
    linexpr1 lexpr = linexpr1(env,linexpr0(std::vector<coeff>(0),coeff(0)));
    var dst = var(std::string("layer"+std::to_string(layer)+"dim"+std::to_string(dim)));
    _elt->assign(*_man,dst, texpr1(lexpr));
  }
  
  std::vector<interval> get_bounds(unsigned layer){
    std::vector<interval> re;
    for(int i = 0; i < _dims.size(); i++)
      re.push_back(get_bound(layer,i));
    return re;
  }


  interval get_bound(unsigned layer, unsigned dim){
    var dst = var(std::string("layer"+std::to_string(layer)+"dim"+std::to_string(dim)));
    return _elt->bound(*_man,dst);
  }


  void join(AbsEle *abs){
    _elt->join(*_man, *(abs->_elt));
  }

  AbsEle* copy(){
    std::vector<var> tdims;
    tdims = _dims;
    manager *tman = new manager(*_man);
    abstract1 *telt = new abstract1(*_man,*_elt);
    *telt = *_elt;
    return (new AbsEle(tman,tdims,telt)); 
  }

  friend std::ostream& operator<< (std::ostream& os, const AbsEle& a){
    std::cout<<*(a._elt)<<std::endl;
  }

 private:
  std::vector<var> _dims;
  abstract1 *_elt;
  manager *_man;
};

#endif
