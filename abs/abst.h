/*********************                                                        */
/*! \file abst.h
 ** Top contributors (to current version):
 **   Jiangchao Liu
 ** This file is part of the AI & SMTReluplex project.
 ** Copyright (c) 2018-2100 by the authors listed in the file AUTHORS
 ** (in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **/

#ifndef __ABST_H
#define __ABST_H


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
#include "apronxx/apxx_linexpr1.hh"


#define AFORP Abst_symV

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

mpq_class get_mpq(double t){
  scalar::order tmp;
  return scalar(t).to_mpq(GMP_RNDN,tmp);
}

// caution:: here the MPQ option should be consistent 
// with the linked library   
coeff get_coeff(double t){
  return coeff(t);
}

var genVar(int layer, int dim){
  std::stringstream ss;
  ss<<"layer"<<layer<<"dim"<<dim;
  var tt(ss.str());
  return tt;
}

linexpr1 vec2expr(environment env, int layer, int ndim, double *coffs,  double consts){
    unsigned size = env.get_vars().size();
    coeff *tcoeff = new coeff[size];		
    unsigned realDim;
    for (int j = 0; j < ndim; j++){
      var cur = genVar(layer, j);
      realDim = env.get_dim(cur);
      tcoeff[realDim] = get_coeff(coffs[j]);
     }
     linexpr1  lexpr1(env,linexpr0(size,tcoeff,get_coeff(consts)));
     delete[] tcoeff;
     return lexpr1;
   }

lincons1 vec2lcon(environment env, int layer,  int ndim, double *coffs, 
		    LogicOP log, double consts){ 
  unsigned size = env.get_vars().size();
    coeff *tcoeff = new coeff[size];
    int sgn = logSign[log];
    unsigned realDim;
    for (int j = 0; j < ndim; j++){
      var cur = genVar(layer, j);
      realDim = env.get_dim(cur);
      tcoeff[realDim] = get_coeff(sgn*coffs[j]);
    }
    linexpr0 t_lexpr0(size,tcoeff,get_coeff(-sgn*(consts)));
    lincons0 t_lcons0(logAp[log],t_lexpr0);
    delete[] tcoeff;
    lincons1 lcon1(env, t_lcons0);
    return lcon1;
  }
  

class Abst{
 public:
  Abst(){
  }
  void add_var(var a){
    throw "not implemented";
  }
  interval get_bound(var a){
    throw "not implemented";
  }
  void rem_var(var a){
    throw "not implemented";
  }
  void join(Abst* a){
    throw "not implemented";
  }
  void meet(Abst* a){
    throw "not implemented";
  }
  void rename(var oldV, var newV){
    throw "not implemented";
  }
  void meet_tcons(tcons1 a){
    throw "not implemented";
  }

  void meet_lcons(lincons1 a){
    throw "not implemented";
  }
  bool is_bottom(){
    throw "not implemented";
  }
  void assign(var dst, linexpr1 aexpr){
    throw "not implemented";
  }
  void powerset(){
    throw "not implemented";
  }
  environment get_environment(){
    throw "not implemented";
  }
  std::vector<var> get_vars(){
    throw "not implemented";
  }
  friend std::ostream& operator<< (std::ostream& os, const Abst& a){
    throw "not implemented";
  }
  ~Abst(){
  }
}; 


class Abst_apron: public Abst{
 public:
  Abst_apron(Abst_apron *a){
    _man = new manager(*(a->_man));
    _elt = new abstract1(*_man, *(a->_elt));
  }
  Abst_apron(Domain aDom){
    environment env = environment({},{});
    switch (aDom){
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
    _elt = new abstract1(*_man,env,top());
  }
  void add_var(var a){
    environment env = _elt->get_environment();
    env = env.add({},{a});
    _elt->change_environment(*_man,env);
  }
  void rem_var(var a){
    environment env =_elt->get_environment();
    env = env.remove({a});
    _elt->change_environment(*_man, env);
  }
  void output_bound(Abst_apron* a){
    environment env = a->_elt->get_environment();
    std::vector<var> vars = env.get_vars();
    std::vector<var>::iterator it;
    for(it = vars.begin(); it != vars.end(); it++){
      std::cout<<*it<<" "<<
	double(a->get_bound(*it).get_inf())<<"  "<<
	double(a->get_bound(*it).get_sup())<<"  "<< std::endl;;
    }
  }
  void join(Abst_apron* a){
    environment env = _elt->get_environment();
    environment env2 = a->_elt->get_environment();
    
    if(cmp(env,env2))
      a->_elt->change_environment(*_man,env);
     
    _elt->join(*_man, *(a->_elt));
  }
  void meet(Abst_apron* a){
    environment env = _elt->get_environment();
    environment env2 = a->_elt->get_environment();
    
    /* std::cout<<"env is"<<env<<std::endl; */
    /* std::cout<<"env2 is"<<env2<<std::endl; */
    /* std::cout<<"elt is"<<std::endl;			 */
    /* output_bound(this); */
    /* std::cout<<"a elt is"<<std::endl; */
    /* output_bound(a); */
    if(cmp(env,env2))
      a->_elt->change_environment(*_man,env);
     
    /* std::cout<<"after evn is chaned a elt is"<<std::endl; */
    /* output_bound(a); */
    _elt->meet(*_man, *(a->_elt));
    /* std::cout<<"join result is"<<std::endl; */
    /* output_bound(this); */
  }
  void rename(var oldV, var newV){
    _elt->rename(*_man,{oldV},{newV});
  } 
  void meet_lcons(lincons1 a){
    _elt->meet(*_man,tcons1_array({a}));
  }
  bool is_bottom(){
   return _elt->is_bottom(*_man);
  }
  void assign(var dst, linexpr1 aexpr){
    //std::cout<<dst<<" "<<aexpr<<std::endl;
    //std::cout<<" before assign "<<std::endl;
    //output_bound(this);
    _elt->assign(*_man, dst, aexpr);
    //std::cout<<" after assign "<<std::endl;
    //output_bound(this);
  }
  interval get_bound(var dst){
    return _elt->bound(*_man, dst);
  }
  environment get_environment(){
    return _elt->get_environment();
  }
  std::vector<var> get_vars(){
    environment env = get_environment();
    return env.get_vars();
  }
  void powerset(){
  }
  ~Abst_apron(){
    delete _man;
    delete _elt;
  }
  friend std::ostream& operator<< (std::ostream& os, const Abst_apron& a){
    os<<*(a._elt)<<std::endl;
    return os;
  }
 private:
  manager *_man;
  abstract1 *_elt;
};

typedef std::map<var, double> Expr;
typedef std::map<var, linexpr1> ExprMap;
typedef std::map<var, Expr> LExprMap;
typedef std::map<var, double> LCstMap;


class Abst_symV: public Abst{
 public:
  Abst_symV(Domain aDom){
    _abst = new Abst_apron(aDom);
  }
  Abst_symV(Abst_symV *a){
    _abst = new Abst_apron(a->_abst);
    _renv = a->_renv;
    _exprMap = a->_exprMap;
    _cstMap = a->_cstMap;
    _symv = a->_symv;
  }
  void powerset(){
  }
  void add_var(var a){
    _abst->add_var(a);
    _renv = _renv.add({},{a});
  }
  void rem_var(var a){
    _renv = _renv.remove({a});
    if(_symv.contains(a)) {
      _symv = _symv.remove({a});
      // _abst->rem_var(a);
      _exprMap.erase(a);
      _cstMap.erase(a);
    }
  }

  void join(Abst_symV* a){
    if (a->is_bottom()) return;
    if (_abst->is_bottom()) {
      _renv = a->_renv;
      _exprMap = a->_exprMap;
      _cstMap = a->_cstMap;
      _abst->join(a->_abst);
      _symv = a->_symv;
      return;
    }
    /* std::vector<var> symv = _symv.get_vars(); */
    /* std::vector<var> asymv = a->_symv.get_vars(); */
    /* std::vector<var>::iterator it; */
    /* for(it = symv.begin(); it!=symv.end(); it++){ */
    /*   if(std::find(asymv.begin(), asymv.end(),*it)==asymv.end()){ */
    /* 	_symv = _symv.remove({*it}); */
    /* 	_exprMap.erase(*it); */
    /* 	_cstMap.erase(*it); */
    /*   } */
    /* }  */
    _symv = a->_symv;
    _abst->join(a->_abst);
  } 
  void meet(Abst_symV* a){
    _abst->meet(a->_abst);
  }

  void rename(var oldV, var newV){
    _abst->rename(oldV, newV);
    /* ExprMap::iterator it = _exprMap.find( oldV ); */
    /* if( it != _exprMap.end() ){ */
    /*   _exprMap.insert(std::pair<var,linexpr1>(newV,it->second)); */
    /*   _exprMap.erase(it); */
    /* } */
    /* for(it = _exprMap.begin(); it != _exprMap.end(); it++){ */
    /*   Expr expr = getExpr(it->second); */
    /*   Expr::iterator tit = expr.find( oldV ); */
    /*   if( tit != expr.end() ) { */
    /* 	expr[newV] = expr[oldV]; */
    /* 	expr.erase(tit); */
    /* 	linexpr1 lexpr(getLinexpr(expr, it->second.get_cst())); */
    /* 	_exprMap.insert(std::pair<var,linexpr1>(it->first, lexpr)); */
    /*   } */
    /* } */
  }
  
  void meet_lcons(lincons1 a){
    lincons1 lcons(a, _renv);
    _abst->meet_lcons(lcons);
  }
  
  void meet_tcons(tcons1 a){
    tcons1 tcons(a, _renv);
    _abst->meet_tcons(tcons);
  }
  
  bool is_bottom(){
    return _abst->is_bottom();
  }

  void assign(var dst, linexpr1 aexpr){
    Expr expr = getExpr(aexpr);
    Expr::iterator bit;
    bool flag = true;
    for(bit = expr.begin(); bit != expr.end(); bit++)
      if(bit->second > 0.000000001 || bit->second <-0.0000000001)
	{flag = false; break;}
    if(flag){
      _abst->assign(dst,aexpr);
      _exprMap.erase(dst);
      _cstMap.erase(dst);
      if (_symv.contains(dst)) {
	_symv = _symv.remove({dst});
      }
      return;
    }
    environment env = _abst->get_environment();
    linexpr1 extExpr(aexpr, env);
    Expr acc_expr = expr;
    double cst = double(extExpr.get_cst().get_scalar());
    Expr::iterator it;
    for(it = expr.begin(); it != expr.end(); it++){
      //here skip the 0 coeff cases
      LExprMap::iterator eit = _exprMap.find(it->first);
      if(eit != _exprMap.end() && _symv.contains(it->first)){
	Expr texpr = eit->second;
	Expr::iterator tit;
	for(tit = texpr.begin(); tit != texpr.end(); tit++){
	  Expr::iterator ttit = acc_expr.find(tit->first);
	  if(ttit != acc_expr.end())
	    acc_expr[tit->first] = acc_expr[tit->first]
	      + it->second * texpr[tit->first];
	  else
	    acc_expr[tit->first] = it->second *texpr[tit->first];
	}
	cst = cst + it->second * _cstMap[it->first] ; 
	acc_expr[it->first] = 0;
      }
    }
    linexpr1 fexpr(getLinexpr(acc_expr, cst));
    _exprMap.insert(std::pair<var, Expr>(dst, acc_expr));
    _cstMap.insert(std::pair<var, double>(dst, cst));
    if (!_symv.contains(dst))  _symv = _symv.add({},{dst});
    _abst->assign(dst, fexpr);
    /* scalar tmp(1); */
    /* lincons1 ori(AP_CONS_EQ,aexpr,tmp); */
    /* _abst->meet_lcons(ori); */
  }
  
  interval get_bound(var dst){
    return _abst->get_bound(dst);
  }
  
  environment get_environment(){
    return _renv;

  }
  std::vector<var> get_vars(){
    return _renv.get_vars();
  }
  ~Abst_symV(){
    delete _abst;
  }
  friend std::ostream& operator<< (std::ostream& os, const Abst_symV& a){
    //LExprMap::iterator it;
    // LExprMap em = a._exprMap;
    //  os<<"the expr map is"<<std::endl;
    //for (it = em.begin(); it != em.end(); it++)
    //  os<<it->first<<" "<<it->second<<std::endl;
    return os<<*(a._abst)<<std::endl;
  }
 private:
  Expr getExpr( linexpr1 lin1 ){
    environment env = lin1.get_environment();
    std::vector<var> vars = env.get_vars();
    std::vector<var>::iterator it;
    Expr expr;
    for(it = vars.begin(); it != vars.end(); it ++){ 
      double value = lin1[*it].get_scalar();
      if(value>0.000000001 || value < -0.00000001) 
	expr[*it] = value;
  }

    return expr;
  }
  linexpr1 getLinexpr(Expr expr, double cst){
    environment env = _abst->get_environment();
    unsigned size = env.intdim() + env.realdim();
    coeff* coffs = new coeff[size];
    Expr::iterator it;
    //std::cout<<"env is \n"<<env<<std::endl;
    for(it = expr.begin(); it != expr.end(); it++){
      // std::cout<<"it is \n"<<it->first<<std::endl; 
      unsigned index = env[it->first];
      coffs[index] = get_coeff(it->second);
    }
    linexpr0 lin0(size, coffs, get_coeff(cst));
    delete[] coffs;
    return linexpr1(env, lin0);
  }
  
  LExprMap _exprMap;
  LCstMap _cstMap;
  environment _renv;
  environment _symv;// It records removable variables
  Abst_apron *_abst;
  
 };


class Abst_ps{
 public:
  Abst_ps(Domain aDom,  unsigned num = 5){
    _num = num;
    _curEleNum = 1;
    AFORP* abst = new AFORP(aDom);
    _abstVec.push_back(abst);
  }
  Abst_ps(Abst_ps *a){
    _num = a->_num;
    _curEleNum = a->_curEleNum;
    AFORP *abst;
    for(unsigned i = 0; i < _curEleNum; i++){
      abst = new AFORP(a->_abstVec[i]);
      _abstVec.push_back(abst);
    }

  }
  void add_var(var a){
    for(unsigned i = 0; i < _curEleNum; i++)
      _abstVec[i]->add_var(a);
  }
  void rem_var(var a){
    for(unsigned i = 0; i < _curEleNum; i++)
      _abstVec[i]->rem_var(a);
  }

  void join(Abst_ps* a){
    for( unsigned j = 0; j < a->_curEleNum;  j++){
      if(!(a->_abstVec[j]->is_bottom()))
	_abstVec[j]->join(a->_abstVec[j]);
    }
  }
  

  void meet(Abst_ps* a){
      for( unsigned j = 0; j < a->_curEleNum;  j++){
	if(!(a->_abstVec[j]->is_bottom()))
	  _abstVec[j]->join(a->_abstVec[j]);
      }
  }
  

  void rename(var oldV, var newV){
    for(unsigned i = 0; i < _curEleNum; i++)
      _abstVec[i]->rename(oldV, newV);
  }
  
  void meet_lcons(lincons1 a){
    for(unsigned i = 0; i < _curEleNum; i++)
      _abstVec[i]->meet_lcons(a);
  }
  
  void meet_tcons(tcons1 a){
    for(unsigned i = 0; i < _curEleNum; i++)
      _abstVec[i]->meet_tcons(a);
  }
  
  bool is_bottom(){
    bool flag = true;
    for(unsigned i = 0; i < _curEleNum; i++)
      if(! _abstVec[i]->is_bottom()) {flag = false; break;}
    return flag;
  }

  void assign(var dst, linexpr1 aexpr){
    for(unsigned i = 0; i < _curEleNum; i++)
      _abstVec[i]->assign(dst, aexpr);
  }
  
  interval get_bound(var dst){
    scalar low = scalar(infty(1));
    scalar high = scalar(infty(-1));
    //for(unsigned i = 0; i < _curEleNum; i++){
    for(unsigned i = 0; i < 1; i++){
      interval tmp = _abstVec[i]->get_bound(dst);
      if( low > tmp.get_inf()) low = tmp.get_inf();
      if( high < tmp.get_sup()) high = tmp.get_sup();
    }						
    return interval(low, high);
  }
  
  environment get_environment(){
    return _abstVec[0]->get_environment();

  }
  std::vector<var> get_vars(){
    return _abstVec[0]->get_environment().get_vars();
  }
  void powerset(){
    environment env = _abstVec[0]->get_environment();
    std::vector<var> vars = env.get_vars();
    assert(_num <= vars.size());
    double *coff = new double[vars.size()];
    unsigned sum = 1;
    for(unsigned i = 0; i < _num; i++){
      std::cout<<"I am running"<<std::endl;
      interval bound = _abstVec[0]->get_bound(vars[i]);
      double med = (double(bound.get_inf()) 
		    + double(bound.get_sup())) / 2;
      sum = 1;
      for(unsigned j = 0; j < i; j++) sum *= 2;
      for(unsigned k = 0; k < sum/2 ; k++){
	AFORP* abst =  new AFORP(_abstVec[k]);
	_abstVec.push_back(abst);
	memset(coff, 0, vars.size() * sizeof(double));
	coff[i] = 1;
	LogicOP log = ABS_SUPEQ;
	lincons1 cons = vec2lcon(env,0,vars.size(), coff, log, med);
	_abstVec[k]->meet_lcons(cons);
	memset(coff, 0, vars.size() * sizeof(double));
	coff[i] = 1;
	log = ABS_LOWEQ;
	cons = vec2lcon(env,0,vars.size(), coff, log, med);
	_abstVec[k+sum/2]->meet_lcons(cons);
      }
    }
    _curEleNum = sum;
    delete[] coff;
  }
  
  ~Abst_ps(){
    for(unsigned i = 0; i < _curEleNum; i++){
      delete _abstVec[i];
    }
    
  }
  friend std::ostream& operator<< (std::ostream& os, const Abst_ps& a){
    os<<"num "<<a._num<<" cur "<<a._curEleNum<<std::endl;
    for(unsigned i = 0; i < a._curEleNum; i++){
      os<<*(a._abstVec[i])<<std::endl;
    }
    return os;
  }
  
 private:
  unsigned _num;
  unsigned _curEleNum;
  std::vector<AFORP*> _abstVec;
  
};





#endif
