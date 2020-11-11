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

#ifndef __LAYER_H
#define __LAYER_H


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "absEle.h"

class Layer{
 public:
  Layer(){
    _layerNum = -1;
  }
  void virtual trans(unsigned layer, AbsEle *ele){
  }
  void set(int layer, int idim, int odim, double **weights, double **bias){
  }
  ~Layer(){
  }

protected:
  int _layerNum;
};

class FCLayer:public  Layer{
 public:
  FCLayer(){
  }
  
  FCLayer(int layer, int idim, int odim, double **weights, double **bias){
    _layerNum = layer;
    _idim = idim;
    _odim = odim;
    _weights = weights;
    _bias = bias;
  } 

  void set(int layer, int idim, int odim, double **weights, double **bias){
    _layerNum = layer;
    _idim = idim;
    _odim = odim;
    _weights = weights;
    _bias = bias;
  }
  
  void RELU(unsigned layer, AbsEle * ele){
    for(int i = 0; i < _odim; i++){
      AbsEle *tmp = ele->copy();
      //  AbsEle *tmp2 = ele->copy();
      ele->assert_positive(layer, i, _odim);
      // tmp2->assert_positive(layer, i, _odim);
      /* interval bound = ele->get_bound(layer, i); */
      //     std::cout<<" layer "<<layer<<" dim "<<i<<" ";
      /* std::cout<< std::setprecision (std::numeric_limits<double>::digits10 + 1) */
      /* 	  << double(bound.get_inf()) <<'\t' */
      /* 	  << double(bound.get_sup()) <<'\t'<<std::endl; */
      
   
      tmp->assert_negative_strict(layer, i, _odim);
      /* bound = tmp->get_bound(layer, i); */
      /* std::cout<<"tmp layer"<<layer<<"dim"<<i<<" "; */
      /* std::cout<< std::setprecision (std::numeric_limits<double>::digits10 + 1) */
      /* 	  << double(bound.get_inf()) <<'\t' */
      /* 	  << double(bound.get_sup()) <<'\t'<<std::endl; */
      
      tmp->assign_zero(layer, i);
      /* bound = tmp->get_bound(layer, i); */
      /* std::cout<<"tmp layer"<<layer<<"dim"<<i<<" "; */
      /* std::cout<< std::setprecision (std::numeric_limits<double>::digits10 + 1) */
      /* 	  << double(bound.get_inf()) <<'\t' */
      /* 	  << double(bound.get_sup()) <<'\t'<<std::endl; */
      
      
      ele->join(tmp);
      //tmp2->join(tmp);
      /* bound = ele->get_bound(layer, i); */
      /* std::cout<<"after join layer"<<layer<<"dim"<<i<<" "; */
      /* std::cout<< std::setprecision (std::numeric_limits<double>::digits10 + 1) */
      /* 	  << double(bound.get_inf()) <<'\t' */
      /* 	  << double(bound.get_sup()) <<'\t'<<std::endl; */
      //delete ele;
      delete tmp;
      // ele = tmp2;
    }
  }

  void linRELU(unsigned layer, AbsEle *ele){
    for(int i = 0; i < _odim; i++){
      ele->linRELU(layer,i);
    }
  }
  
  
  
  void trans(unsigned layer, AbsEle *ele){
   assert(ele);
   //   std::cout<<"layer "<<layer<<std::endl;
   ele->matrix_assign(layer,_idim,_odim,_weights,_bias);
   std::vector<interval> bounds = ele->get_bounds(layer, _odim);
   std::ofstream ofre("bounds.txt",std::ios::app);
   //   ofre << "*" << _layerNum <<std::endl;
   for(unsigned i = 0; i < bounds.size();i++)
       ofre //<< std::scientific
            << std::setprecision (std::numeric_limits<double>::digits10 + 1)
            << double(bounds[i].get_inf()) <<'\t'
            << double(bounds[i].get_sup()) <<'\t'<<std::endl;
   ofre<<std::endl;
   //std::cout<<*ele<<std::endl;
   RELU(layer, ele);
   bounds = ele->get_bounds(layer, _odim);
   /* ofre << "***" << _layerNum <<std::endl; */
   /* for(unsigned i = 0; i < bounds.size();i++) */
   /*     ofre //<< std::scientific */
   /*          << std::setprecision (std::numeric_limits<double>::digits10 + 1) */
   /*          << double(bounds[i].get_inf()) <<'\t' */
   /*          << double(bounds[i].get_sup()) <<'\t'<<std::endl; */
   /* ofre<<std::endl; */
   ofre.close();
  }

 private:
  int _idim;
  int _odim;
  double **_weights;
  double **_bias;
};
#endif
