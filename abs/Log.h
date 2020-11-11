#ifndef __LOG_H
#define __LOG_H 

#include<iostream>
#include<string>
void Log_info(const char* str){
  std::cout<<str<<std::endl<<std::flush;
}
#endif
