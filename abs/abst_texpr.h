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
#include "apronxx/apxx_texpr1.hh"
#include "apronxx/apxx_tcons1.hh"

using namespace apron;
enum Domain {BOX, OCT, POLKA, TOPO};

class Abst{
 public:
  Abst(){
    std::cout<<"construction upper layer"<<std::endl;
  }
  void add_var(var a){
    std::cout<<"I should not be called"<<std::endl;
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
  void rename(var oldV, var newV){
    throw "not implemented";
  }
  void meet_tcons(tcons1 a){
    throw "not implemented";
  }
  bool is_bottom(){
    throw "not implemented";
  }
  void assign(var dst, texpr1 aexpr){
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
    std::cout<<"construction lower layer"<<std::endl;
    _man = new manager(*(a->_man));
    _elt = new abstract1(*(a->_man), *(a->_elt));
  }
  Abst_apron(Domain aDom){
    std::cout<<"construction lower layer"<<std::endl;
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
    _elt->change_environment(*_man,env);
  }
  void join(Abst_apron* a){
    _elt->join(*_man, *(a->_elt));
  }
  void rename(var oldV, var newV){
    _elt->rename(*_man,{oldV},{newV});
  } 
  void meet_tcons(tcons1 a){
    _elt->meet(*_man,tcons1_array({a}));
  }
  bool is_bottom(){
   return _elt->is_bottom(*_man);
  }
  void assign(var dst, texpr1 aexpr){
    _elt->assign(*_man, dst, aexpr);
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
  ~Abst_apron(){
    delete _man;
    delete _elt;
  }
  friend std::ostream& operator<< (std::ostream& os, const Abst_apron& a){
    std::cout<<*(a._elt)<<std::endl;
    return os;
  }
 private:
  manager *_man;
  abstract1 *_elt;
};

/* typedef std::map<var, coeff> Texpr; */
/* typedef std::map<var, Texpr> Tsyms; */


/* class Abst_symV: public Abst{ */
/*  public: */
/*   Abst_symV(Domain aDom){ */
/*     _abst = new Abst_apron(aDom); */
/*   } */
/*   void add_var(var a){ */
/*     _abst->add_var(a); */
/*     _renv.push_back(a); */
/*   } */
/*   void rem_var(var a){ */
/*     _renv.erase(a); */
/*     clean_var(a); */
/*   } */

/*   void join(Abst_symV* a){ */
/*     _abst->join(*a); */
/*     Tsyms::iterator it; */
/*     for(it = _symbols.begin(); it != _symbols.end(); it++){ */
/*       Tsyms::iterator tit = a->_symbols.find(it->first); */
/*       if(tit == a->_symbols.end()) _symbols.erase(it); */
/*       else if(it->second != tit->second) _symbols.erase(it); */
/*     } */
    
/*   } */
/*   void rename(var oldV, var newV){ */
/*     _abst->rename(oldV, newV); */
/*     Tsyms::iterator it = _symbols.find( oldV ); */
/*     if( it != _symbols.end() ){ */
/*       _symbols[newV] = it->second; */
/*       _symbols.erase(it); */
/*     } */
/*     for(it = _symbols.begin(); it != _symbols.end(); it++){ */
/*       Texpr::iterator tit = it->second.find( oldV ); */
/*       if( tit != it->second.end() ) { */
/* 	it->second[newV] = tit->second; */
/* 	it->second.erase( tit ); */
/*       } */
/*     } */
/*   }  */
  
/*   void meet_tcons(tcons1 a){ */
/*     _abst->meet_tcons(a); */
/*   } */
  
/*   bool is_bottom(){ */
/*     _abst->is_bottom(); */
/*   } */

/*   void assign(var dst, texpr1 aexpr){ */
/*     if (aexpr.is_interval_cst()){ */
/*       _abst->assign(dst, aexpr); */
/*       _symbols.erase(dst); */
/*       return; */
/*     } */
/*     Texpr te = getTexpr( const_iterator( aexpr ) ); */
/*     Texpr::iterator it; */
/*     for( it = te.begin(); it != te.end(); it++){ */
/*       Tsyms::iterator tit = _symbols.find(it->first); */
/*       if (tit != _sysmbols.end() && tit->first != dst){ */
/* 	Texpr src = tit->second; */
/* 	Texpr::iterator eit; */
/* 	for( eit = src.begin(); eit != src.end(); eit++){ */
/* 	  Texpr::iterator etit = te.find(eit->first); */
/* 	  if(etit != te.end()) te[eit->first] = te[eit->first] + eit->second; */
/* 	  else te[eit->first] = eit->second; */
/* 	} */
/*       } */
/*     } */
/*     symbols[dst]=te; */
/*     _abst->assign(dst, get_texpr1(te)); */
/*   } */

/*   interval get_bound(var dst){ */
/*     _abst->get_bound(dst); */
/*   } */

/*   environment get_environment(){ */
/*     return environment({},_renv); */

/*   } */
/*   std::vector<var> get_vars(){ */
/*   } */
/*   ~Abst_apron(){ */
/*   } */
/*   friend std::ostream& operator<< (std::ostream& os, const Abst_symV& a){ */
/*   } */
/*  private: */
/*   Texpr getTexpr( const_iterator it ){ */
/*     siwtch() */
/*   } */
/*   void clean_var(var a){ */
/*     Tsyms::iterator it; */
/*     bool flag = true; */
/*     for(it = _symbols.begin(); it != _symbols.end(); it++){ */
/*       Texpr aex = it->second; */
/*       Texpr::iterator tit = aex.find( a ); */
/*       if(tit != aex.end()){ */
/* 	flag = false; */
/* 	break; */
/*       } */
/*     } */
/*     if(flag) _abst->rem_var(a);   */
/*   } */
/*   void cleanup(){ */
/*     std::vector<var> vars = _abst->get_vars(); */
/*     std::vactor<var>::iterator it; */
/*     for (it = vars.begin(); it != vars.end(); it++){ */
/*       std::vactor<var>::iterator tit = _renv.find(it); */
/*       if(tit == _renv.end()) clean_var(*it); */
/*     } */
/*   } */
  
/*   Tsyms _symbols; */
/*   std::vector<var> _renv; */
/*   Abst_apron *_abst; */
  
/*  }; */



typedef std::map<var, texpr1> Tsyms;


class Abst_symV: public Abst{
 public:
  
  Abst_symV(Abst_symV* as){
    _abst = new Abst_apron(as->_abst);
    _renv = environment(as->_renv);
    _symbols = as->_symbols;
  }
  
  Abst_symV(Domain aDom){
    _abst = new Abst_apron(aDom);
    _renv = environment();
  }
  void add_var(var a){
    _abst->add_var(a);
    _renv = _renv.add({},{a});
  }
  void rem_var(var a){
    _renv = _renv.remove({a});
    clean_var(a);
  }
  
  bool texpr1Eq(texpr1::const_iterator a, texpr1::const_iterator b){
    return true;
    if ( a.get_discr() != b.get_discr() ) return false;
    switch(a.get_discr()){
    case AP_TEXPR_CST:
      return a.get_coeff() == b.get_coeff();
    case AP_TEXPR_DIM:
      return a.get_var() == b.get_var();
    case AP_TEXPR_NODE:
      switch(a.get_op()){
      case AP_TEXPR_NEG:
      case AP_TEXPR_CAST:
      case AP_TEXPR_SQRT:
  	return a.get_op()==b.get_op() && texpr1Eq(a.left(), b.left());
      default:
  	return a.get_op()==b.get_op()
  	  && texpr1Eq(a.left(), b.left())
  	  && texpr1Eq(a.right(), b.right());
      }
    default:
      throw "unexpected case";
    }
  }
  
  void replace(texpr1::iterator& it, var dst, texpr1::builder& x){
    switch(it.get_discr()){
    case AP_TEXPR_CST:
      return;
    case AP_TEXPR_DIM:
      if( it.get_var() == dst){
  	it = x;
      }
      return;
    case AP_TEXPR_NODE:
      switch(it.get_op()){
      case AP_TEXPR_NEG:
      case AP_TEXPR_CAST:
      case AP_TEXPR_SQRT:{
	texpr1::iterator tit(it.left());
  	replace(tit, dst, x);
  	return;
      }
      default:{
	texpr1::iterator tit0(it.left());
	texpr1::iterator tit1(it.right());
  	replace(tit0, dst, x);
  	replace(tit1, dst, x);
  	return;
      }
      }
    default:
      throw "unexpected case";
    }
  }


  void join(Abst_symV* a){
    _abst->join(a->_abst);
    Tsyms::iterator it;
    for(it = _symbols.begin(); it != _symbols.end(); it++){
      Tsyms::iterator tit = a->_symbols.find(it->first);
      if(tit == a->_symbols.end()) _symbols.erase(it);
      else if(!texpr1Eq(texpr1::const_iterator((it->second)),
  			texpr1::const_iterator((tit->second))))
  	_symbols.erase(it);
    }
    cleanup();
  }
  void rename(var oldV, var newV){
    _abst->rename(oldV, newV);
    Tsyms::iterator it = _symbols.find( oldV );
    if( it != _symbols.end() ){
      _symbols.insert(std::pair<var,texpr1>(newV,it->second));
      _symbols.erase(it);
    }
    for(it = _symbols.begin(); it != _symbols.end(); it++){
      texpr1::iterator tit(it->second);
      texpr1::builder b(texpr1::builder(_abst->get_environment(),newV));
      replace(tit, oldV, b);
    }
  }
  
  void meet_tcons(tcons1 a){
    _abst->meet_tcons(a);
  }
  
  bool is_bottom(){
    return _abst->is_bottom();
  }

  void assign(var dst, texpr1 aexpr){
    if (aexpr.is_interval_cst()){
      std::cout<<"this should not be called"<<std::endl;
      _abst->assign(dst, aexpr);
      _symbols.erase(dst);
      return;
    }
    Tsyms::iterator it;
    for(it = _symbols.begin(); it != _symbols.end(); it++){
      std::cout<<"aexpr is "<<aexpr<<std::endl;
      std::cout<<"first is "<<it->first<<std::endl;
      if (aexpr.has_var(it->first)) {
  	texpr1 ex(it->second);
  	texpr1::builder b(ex);
	texpr1::iterator tit(aexpr);
	std::cout<<"replace "<<it->first<<" to "<<b<<std::endl;
  	replace(tit, dst, b);
      }
    }
    _symbols.insert(std::pair<var,texpr1>(dst,aexpr));
    std::cout<<"attention"<<std::endl;
    std::cout<<dst<<"="<<aexpr<<std::endl;
    _abst->assign(dst, aexpr);
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
  }
  friend std::ostream& operator<< (std::ostream& os, const Abst_symV& a){
    std::cout<<*(a._abst)<<std::endl;
    return os;
  }

 private:
  
  void clean_var(var a){
    Tsyms::iterator it;
    bool flag = true;
    for(it = _symbols.begin(); it != _symbols.end(); it++){
      texpr1 aex(it->second);
      if( aex.has_var(a)){
  	flag = false;
  	break;
      }
    }
    if(flag) {
      _abst->rem_var(a);
      _symbols.erase(a);
    }
  }
  void cleanup(){
    std::vector<var> vars = _abst->get_vars();
    std::vector<var>::iterator it;
    for (it = vars.begin(); it != vars.end(); it++){
      if(!_renv.contains(*it)) clean_var(*it);
    }
  }
  
  Tsyms _symbols;
  environment _renv;
  Abst_apron *_abst;
  
 };



#endif
