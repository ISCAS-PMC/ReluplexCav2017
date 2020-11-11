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

#ifndef __ABSEXE_H
#define __ABSEXE_H


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
#include <boost/date_time.hpp>

int maxRelu = 0;
int inActive = 0;

enum ROBUST {AI, PLANET};
using namespace apron;

void robust_rlv(std::string img, ROBUST trobust, 
		double delta, std::string network, Domain dom,
        string bound_file_name, string summary_file_name
){

    cout << "img file name     \t: \t" << img << endl;
    cout << "rlv net file name \t: \t" << network << endl;
    cout << "bound file name   \t: \t" << bound_file_name << endl;
    cout << "summary file name \t: \t" << summary_file_name << endl;

    cout << "delta             \t= \t" << delta << endl;
    cout << "domain            \t= \t" << (dom == BOX ? "BOX" : "TOPO") << endl;
    cout << "robustness        \t= \t" << (trobust == PLANET ? "PLANET" : "AI2") << endl;

  	std::ofstream ofre(bound_file_name,std::ios::app);
  	std::ifstream inPic(img);
  	std::ifstream rlvNet(network);
  	std::string currentLine;
  	ABST* abst = new ABST(dom);
  	int layer;
  	int input = 0;
	unsigned curNum = 0;
	unsigned nodeNum;
	bool outLayer = false;
	std::vector<interval> boundFinal;
	bool Verified = true;
	auto start_time = boost::posix_time::microsec_clock::local_time();


	while(std::getline(rlvNet, currentLine)){
	  while ((currentLine.size()>0) && (currentLine[0]==' ')) currentLine = currentLine.substr(1,std::string::npos);
	  while ((currentLine.size()>0) && (currentLine[currentLine.size()-1]==' ')) currentLine = currentLine.substr(0,currentLine.size()-1);
         if (currentLine.size()>0) {
	   /* std::cout<<currentLine<<std::endl; */
	   /* std::cout<<"apron assign "<<apron_assign_time<<std::endl; */
	   /* std::cout<<"apron cons "<<apron_cons_time<<std::endl; */
	   /* std::cout<<"apron join "<<apron_join_time<<std::endl; */
	  
	   /* std::cout<<"symv assign "<<symv_assign_time<<std::endl; */
	   /* std::cout<<"symv cons "<<symv_cons_time<<std::endl; */
	   /* std::cout<<"symv join "<<symv_join_time<<std::endl; */
	   // Parse this line
	   std::istringstream thisLineStream(currentLine);
	   std::string linePrefix;
	   thisLineStream >> linePrefix;
	    if (linePrefix == "#"){
	      std::string slayer;
	      unsigned layer, num;
	      string ty, tn;
	      thisLineStream >> slayer >> layer >> nodeNum >> ty >> tn;
	      std::cout<<slayer<<layer<<std::endl;	
	      ofre<<"* "<<layer<<std::endl;
	      if(ty == "Linear" && tn == "Accuracy") outLayer = true;
	    } else if (linePrefix=="Input") {
	      // Input nodes
	      input++;
	      std::string nodeName;
	      double value, up, low;
	      inPic>>value;
	      thisLineStream >> nodeName;
	      var dstVar(nodeName);
	      abst->add_var(dstVar);
	      value = value * 0.00390625;
	      switch (trobust){
	      case AI:
		 low = value;
		 if(value < 1 - delta)
		   up = value;
		 else up = 1;
		 break;
	      case PLANET:
		if (input <= 28 *3 || input >= 28 * 26 
		    || input - input/28 * 28 <= 3 || input - input/28 * 28 >= 26){
		  up = value;
		  low = value;
		}
		else{
		  up =  value + delta < 1 ?value+delta:  1;
		  low =  value - delta < 0? 0 : value - delta;
		}
		break;
	      default:
		throw "unknow robust type";
	      }
	      environment env = abst->get_environment();
	      unsigned  dim = env.get_dim(dstVar);
	      unsigned size = env.get_vars().size();
	      coeff *tcoeff = new coeff[size];
	      tcoeff[dim] = get_coeff(1);
	      linexpr0 t_lexpr0(size,tcoeff,get_coeff(-low));
	      lincons0 t_lcons0(AP_CONS_SUPEQ,t_lexpr0);
	      lincons1 lcons1(env, t_lcons0);
	      abst->meet_lcons(lcons1);
	      tcoeff[dim] = get_coeff(-1);
	      t_lexpr0 = linexpr0(size,tcoeff,get_coeff(up));
	      t_lcons0 = lincons0(AP_CONS_SUPEQ,t_lexpr0);
	      lcons1 = lincons1(env, t_lcons0);
	      abst->meet_lcons(lcons1);	
	      delete[] tcoeff;
            } else if (linePrefix=="MaxPool") {
	      std::string nodeName;
	      thisLineStream >> nodeName;
	      var dstVar = var(nodeName);
	      abst->add_var(dstVar);
	      std::vector<var> pvars;
	      while (!thisLineStream.eof()) {
		std::string src;
		thisLineStream >> src;
		pvars.push_back(var(src));
	      }
	      std::vector<var>::iterator it;
	      std::vector<var>::iterator it2;
	      environment env = abst->get_environment();
	      unsigned size = env.get_vars().size();
	      coeff *tcoeff = new coeff[size];
	      coeff *tcoeffA = new coeff[size];
	      ABST* ori = abst;
	      for(it = pvars.begin(); it != pvars.end(); it++){
		ABST* tmp = new ABST(ori);
		unsigned dim1 = env[*it];
		tcoeff[dim1] = get_coeff(1);
		interval boundO = ori->get_bound(*it);
		double infO = double(boundO.get_inf());
		double supO = double(boundO.get_sup());
		for(it2 = pvars.begin(); it2 != pvars.end(); it2++){
		  if(it == it2) continue;
		  interval boundi = ori->get_bound(*it2);
		  double infi = double(boundi.get_inf());
		  double supi = double(boundi.get_sup());
		  if(infO >= supi) continue;
		  unsigned dim2 = env[*it2];
		  tcoeff[dim2] = get_coeff(-1);
		  linexpr0 t_lexpr0(size,tcoeff,get_coeff(0));
		  lincons0 t_lcons0(AP_CONS_SUPEQ,t_lexpr0);
		  lincons1 lcon1(env, t_lcons0);
		  tmp->meet_lcons(lcon1);
		  tcoeff[dim2] = get_coeff(0);
		  if(tmp->is_bottom()) break;
		}
		tcoeff[dim1] = get_coeff(0);
		if(!(tmp->is_bottom())){
		  tcoeffA[dim1] = get_coeff(1);
		  linexpr0 a_lexpr0(size, tcoeffA,get_coeff(0));
		  linexpr1 a_lexpr1(env,a_lexpr0);
		  tmp->assign(dstVar, a_lexpr1);
		}
		if (it != pvars.begin()) {
		  abst->join(tmp);
		  delete tmp;
		}
		else abst = tmp;
		tcoeffA[dim1] = get_coeff(0);
	      }
	      delete[] tcoeff;
	      delete[] tcoeffA;
	      delete ori;
	      ofre<<nodeName<<" ";
	      interval fbound = abst->get_bound(dstVar);
	      ofre << std::setprecision (std::numeric_limits<double>::digits10 + 1)
		   << double(fbound.get_inf()) <<'\t'
		   << double(fbound.get_sup()) <<'\t'<<std::endl;
	    } else if ( linePrefix=="Linear") {
	      std::string nodeName;
	      thisLineStream >> nodeName;
	      var dstVar(nodeName);
	      abst->add_var(dstVar);
	      double cst;
	      thisLineStream >> cst;
	      double weight;
	      environment env = abst->get_environment();
	      unsigned size = env.get_vars().size();
	      coeff *tcoeff = new coeff[size];
	      while (!thisLineStream.eof()) {
                    thisLineStream >> weight;
		    std::string src;
                    thisLineStream >> src;
                    tcoeff[env[var(src)]] = get_coeff(weight);
	      }
	      linexpr0 a_lexpr0(size, tcoeff,get_coeff(cst));
	      linexpr1 a_lexpr1(env,a_lexpr0);
	      abst->assign(dstVar,a_lexpr1);
	      delete[] tcoeff;
	      ofre<<nodeName<<" ";
	      interval bound = abst->get_bound(dstVar);
	      ofre << std::setprecision (std::numeric_limits<double>::digits10 + 1)
		   << double(bound.get_inf()) <<'\t'
		   << double(bound.get_sup()) <<'\t'<<std::endl;
	      if( outLayer ){
		curNum ++;
		boundFinal.push_back(bound);
		if(curNum == nodeNum){
		  std::vector<interval>::iterator it = boundFinal.begin();
		  double lb = double(it->get_inf());
		  double ub =  double(it->get_sup());
		  for(it = it + 1; it != boundFinal.end(); it++){
		    if( Verified && lb <  double(it->get_sup())) Verified = false; 
		    if(! Verified){
		      if (ub < double(it->get_inf())){
			ub = double(it->get_sup());
			lb = double(it->get_inf());
			Verified = true;
		      }
		      else{
			ub = ub >  double(it->get_sup())? ub:double(it->get_sup()); 
		      }
		    }
		  }
		}
	      }
	      
	    } else if (linePrefix=="ReLU" ) {
	      maxRelu ++;
	      std::string nodeName;
	      thisLineStream >> nodeName;
	      var dstVar(nodeName);
	      abst->add_var(dstVar);
	      double cst;
	      thisLineStream >> cst;
	      double weight;
	      environment env = abst->get_environment();
	      unsigned size = env.get_vars().size();
	      coeff *tcoeff = new coeff[size];
	      while (!thisLineStream.eof()) {
                    thisLineStream >> weight;
		    std::string src;
                    thisLineStream >> src;
                    tcoeff[env[var(src)]] = get_coeff(weight);
	      }
	      linexpr0 a_lexpr0(size, tcoeff,get_coeff(cst));
	      linexpr1 a_lexpr1(env,a_lexpr0);
	      
	      abst->assign(dstVar,a_lexpr1);
	      delete[] tcoeff;
	      /* std::string nodeName; */
	      /* thisLineStream >> nodeName; */
	      /* var dstVar(nodeName); */
	      /* double cst; */
	      /* std::string src; */
	      /* thisLineStream >> cst >> cst >> src; */
	      /* var srcVar(src); */
	      /* abst->rename(srcVar,dstVar); */
	      /* environment env = abst->get_environment(); */
	      /* unsigned size = env.get_vars().size(); */
	      interval bound = abst->get_bound(dstVar);
	      double inf = double(bound.get_inf());
	      double sup =  double(bound.get_sup()); 
	      ofre<<nodeName<<" ";
	      ofre << std::setprecision (std::numeric_limits<double>::digits10 + 1)
		   << inf <<'\t'
		   << sup <<'\t'<<std::endl;	    
	      if(inf >= 0) { inActive ++; continue;}
	      if(sup < 0){
		inActive ++;
		coeff *coff = new coeff[size];
		linexpr0 t_lexpr0 = linexpr0(size,coff,get_coeff(0));
		abst->assign(dstVar, linexpr1(env,t_lexpr0));
		delete[] coff;
		continue;
	      }
	      ABST* tmp = new ABST(abst);
	      coeff *coff = new coeff[size];
	      coff[env[dstVar]] = get_coeff(1);
	      linexpr0 t_lexpr0(size,coff,get_coeff(0));
	      lincons0 t_lcons0(AP_CONS_SUPEQ,t_lexpr0);
	      abst->meet_lcons(lincons1(env, t_lcons0));
	      coff[env[dstVar]] = get_coeff(-1);
	      t_lexpr0 = linexpr0(size,coff,get_coeff(0));
	      t_lcons0 = lincons0(AP_CONS_SUPEQ,t_lexpr0);
	      tmp->meet_lcons(lincons1(env, t_lcons0));
	      if(! tmp->is_bottom()){
		coff[env[dstVar]] = get_coeff(0);
		t_lexpr0 = linexpr0(size,coff,get_coeff(0));
		tmp->assign(dstVar, linexpr1(env,t_lexpr0));
		abst->join(tmp);
	      }
	      delete tmp;
	      delete[] coff;
	    }
	    else if (linePrefix=="Assert") {
	      std::cout<<"assert encountered, but we do nothing"<<std::endl;
	    } else {
                std::ostringstream err;
		std::cout<<"Error: Did not understand line prefix " << linePrefix;
                err << "Error: Did not understand line prefix " << linePrefix;
                throw err.str();
        }

	 }
  }
  auto end_time = boost::posix_time::microsec_clock::local_time();
  auto time_d = end_time - start_time;
  std::cout << time_d << std::endl;


	std::ofstream summary(summary_file_name,std::ios::app);
	summary
			<< network << ",\t"
			<< img << ",\t"
			<< delta << ",\t"
			<< (trobust == PLANET ? "PLANET" : "AI2") << ",\t"
	        << (dom == BOX ? "BOX" : "TOPO") << ",\t"

#ifdef NUM_D
			<< "NUM_D,\t"
#endif

#ifdef NUM_MPQ
	        << "NUM_MPQ,\t"
#endif

#ifdef ABST
	        << ABST::type_name() << ",\t"
#endif
            << time_d << ",\t"
			<< time_d.total_milliseconds() << ",\t"
            << end_time << ",\t"<<maxRelu << ",\t"<<inActive<< ",\t"<<Verified
			<< endl;
}
#endif//__ABSEXE_H
