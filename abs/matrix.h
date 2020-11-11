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
#ifndef __MATRIX_H
#define __MATRIX_H

#include "Log.h"
#include "absEle.h"
#include <iostream>
#include <vector>

using namespace std;

class Matrix{
public:
Matrix(int ndim){			       
_ndim = ndim;

}
void generate_matrix(int dim, double value){
  double *coff = new double[_ndim];
  memset(coff,0,_ndim*sizeof(double));
  coff[dim] = 1;
  _coffs.push_back(coff);
  double *bia = new double[1];
  *bia = value;
  _bias.push_back(bia);
}
 void set_uperbound(int dim, double value){
   generate_matrix(dim,value);
   _lop.push_back(ABS_LOWEQ);
 }
 void set_lowerbound(int dim, double value){
   generate_matrix(dim,value);
   _lop.push_back(ABS_SUPEQ);
 }
 double** get_coff(){
   return _coffs.data();
 }
 LogicOP* get_lop(){
   return _lop.data();
 }
 double** get_bias(){
   return _bias.data();
 }
 int get_num(){
   return _coffs.size();
}
 ~Matrix(){
   for(int i=0;i<_ndim;i++){
     delete[] _coffs[i];
     delete[] _bias[i];
   }
 }


 
 private:
 int _ndim;
 vector<double*> _coffs;
 vector<LogicOP> _lop;
 vector<double*> _bias;
};

#endif
