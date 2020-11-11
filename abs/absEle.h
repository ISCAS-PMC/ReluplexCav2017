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
#include "apronxx/apxx_interval.hh"
#include "apronxx/apxx_oct.hh"
#include "apronxx/apxx_polka.hh"
#include "apronxx/apxx_t1p.hh"
#include "apronxx/apxx_texpr1.hh"
#include "apronxx/apxx_tcons1.hh"
#include "abst.h"
#include "ap_scalar.h"
#include "gmpxx.h"


static int unsplit = 0;
static int split = 0;
using namespace apron;


class AbsEle {

 public:
  
  AbsEle(ABST* abst){
    _abst = abst;
  }
  
  AbsEle(int ndim, std::vector<std::string> name, Domain dom, int num, double** coffs, LogicOP *lop, double **consts){
      
    _abst = new ABST(dom);
     for(int i = 0; i < ndim; i++) 
      _abst->add_var(genVar(0,i));
     for(int i = 0; i < num; i++)
       _abst->meet_lcons(vec2lcon(_abst->get_environment(),
				  0, ndim,coffs[i],lop[i],*(consts[i])));
     _abst->powerset();
     std::cout<<"after initinal\n"<<*_abst<<std::endl;
  }

  ~AbsEle(){
    delete _abst;
  }
  
  void matrix_assign(unsigned layer, int idim, int odim, double **weights, double **bias){
    for(int i = 0; i < odim; i++){
      var dst = genVar(layer,i);
      _abst->add_var(dst);
      linexpr1 expr = vec2expr(_abst->get_environment(),
			       layer-1, idim, weights[i], *bias[i]);
      _abst->assign(dst, expr);
    } 
    for(int i = 0; i < idim; i++)
      _abst->rem_var(genVar(layer - 1, i));
  }
  
  void linRELU(unsigned layer, unsigned dim){
    var cvar = genVar(layer, dim);
    interval bound = _abst->get_bound(cvar);
    if(double(bound.get_sup()) <= 0 ){
      this->assign_zero(layer, dim);
      return;
    }
    if(double(bound.get_inf())>=0) return;
    var tvar(std::string("tmp"));
    _abst->add_var(tvar);
    environment env = _abst->get_environment();
    double  lh = double(bound.get_sup())*double(bound.get_inf());
    double hMinusl = double(bound.get_sup())-double(bound.get_inf());
    double h = double(bound.get_sup());
    unsigned size = env.get_vars().size();
    coeff* tceff = new coeff[size];
    for(unsigned i = 0; i < size; i++) tceff[i] = get_coeff(0);
    tceff[env[cvar]] = get_coeff(h);
    tceff[env[tvar]] = get_coeff(-hMinusl);
    linexpr0 t_lexpr0(size, tceff, get_coeff(-lh));
    lincons0 t_lcons0(AP_CONS_SUPEQ, t_lexpr0);
    lincons1 lcons1(env, t_lcons0);
    _abst->meet_lcons(lcons1);
    
    for(unsigned i = 0; i < size; i++) tceff[i] = get_coeff(0);
    tceff[env[cvar]] = get_coeff(-1);
    tceff[env[tvar]] = get_coeff(1);
    t_lexpr0 = linexpr0(size, tceff, get_coeff(0));
    t_lcons0 = lincons0(AP_CONS_SUPEQ, t_lexpr0);
    lcons1 = lincons1 (env, t_lcons0);
    _abst->meet_lcons(lcons1);
    
    for(unsigned i = 0; i < size; i++) tceff[i] = get_coeff(0);
    tceff[env[tvar]] = get_coeff(1);
    t_lexpr0 = linexpr0(size, tceff, get_coeff(0));
    t_lcons0 = lincons0(AP_CONS_SUPEQ, t_lexpr0);
    lcons1 = lincons1 (env, t_lcons0);
    _abst->meet_lcons(lcons1);
    delete[] tceff;
    _abst->rem_var(cvar);
    _abst->rename(tvar, cvar);
  }


  interval get_bound(unsigned layer, unsigned dim){
    var dst = genVar(layer, dim);
    return _abst->get_bound(dst);
  }

  void assert_positive(int layer, int dim, int odim){
    double *coff = new double[odim];
    memset(coff, 0, odim * sizeof(double));
    coff[dim] = 1;
    LogicOP log = ABS_SUP;
    lincons1 cons = vec2lcon(_abst->get_environment(),layer, odim, coff, log, 0);
    _abst->meet_lcons(cons);
    if (_abst->is_bottom()) unsplit+=2;
    delete[] coff;
  }
  void assert_negative_strict(int layer, int dim, int odim){
    double *coff = new double[odim];
    memset(coff, 0, odim * sizeof(double));
    coff[dim] = 1;
    LogicOP log = ABS_LOWEQ;
    lincons1 cons = vec2lcon(_abst->get_environment(),layer,odim, coff, log, 0);
    _abst->meet_lcons(cons);
    if (_abst->is_bottom()) unsplit++;
    delete[] coff;
  }
  void assign_zero(int layer, int dim){
    if(_abst->is_bottom()) return;
    environment env = _abst->get_environment();
    linexpr1 lexpr = linexpr1(env,linexpr0(std::vector<coeff>(0),get_coeff(0)));
    var dst = genVar(layer, dim);
    _abst->assign(dst, lexpr);
  }
  
  std::vector<interval> get_bounds(unsigned layer, unsigned ndim){
    std::vector<interval> re;
    for(unsigned i = 0; i < ndim ; i++){
      re.push_back(get_bound(layer,i));
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
    return std::cout<<*(a._abst)<<std::endl;
  }

 private:
  ABST* _abst;
};

#endif
