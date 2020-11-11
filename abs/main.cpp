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
#include "Log.h"
#include <string.h>
#include "dnn.h"
#include "absEle.h"
#include "absEngine.h"
#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "matrix.h"
#include "absExe.h"
#include <sys/resource.h>
#include <signal.h>
#include <fstream>
#include <boost/date_time.hpp>


std::string inFileName = "";
std::string summaryFileName = "";
std::string ans = "";

void got_signal( int )
{
    printf( "Got signal\n" );

    if(summaryFileName != "") {

//        auto current = boost::posix_time::microsec_clock::local_time();
//        auto totalTime = current - problem.startTime;
//
//        long long initialTime = problem.initialWorkTime.total_milliseconds();
//        long long solveTime   = problem.solveTime.total_milliseconds();
//        long long total = totalTime.total_milliseconds() ;
//
//
//        std::ofstream summary(summaryFileName,std::ios::app);
//        auto p = inFileName.find_last_of('/');
//        summary << inFileName.substr(p+1,inFileName.length()) << " ,"
//                << "TIMEOUT , "
//                << totalTime << " ,"
//                << total << " ,"
//                << initialTime << " ,"
//                << solveTime << " ,"
//                << problem.cacheSize()<< std::endl;
//        summary.close();
    }


    exit(0);
}


int acas(){
  std::string path("/home/jiangchao/AIProjects/reluplex/ReluplexCav2017/nnet/ACASXU_run2a_1_2_batch_2000.nnet");
  //std::string path("/home/jiangchao/AIProjects/AIsmt/ReluplexCav2017/nnet/test.nnet");
  Log_info("Build an Acas network");
  Dnn *adnn = new Dnn(Acas, path);
  Log_info("Get the dims in the input layer");
  int inputDim = adnn->get_inputdim();
  std::cout<<"dim is"<<inputDim<<std::endl;
  Log_info("Build an engine based on abstract interpretation with an Acas network");
  AbsEngine *aeng = new AbsEngine(adnn);
  std::vector<std::string> name(inputDim);
  int i;
  for(i = 0; i < inputDim; i++){
    name[i] =  "dim" + std::to_string(i);
  }
  Log_info("Introduce contriants on inputs");
  double *mins = adnn->get_mins();
  double *maxes = adnn->get_maxes();
  Matrix amatrix(inputDim);
  // for(int i = 0; i < inputDim; i++){
  //   amatrix.set_uperbound(i,maxes[i]);
  //    amatrix.set_lowerbound(i,mins[i]);   
  // }
  // amatrix.set_lowerbound(0,55947);
  // amatrix.set_lowerbound(3,1145);
  // amatrix.set_uperbound(4,60);

  // amatrix.set_lowerbound(0,-0.2000);
  // amatrix.set_uperbound(0,0.00000);
  // amatrix.set_lowerbound(1,0.027465);
  // amatrix.set_uperbound(1,0.227465);
  // amatrix.set_lowerbound(2,-0.200703);
  // amatrix.set_uperbound(2,-0.000703);
  // amatrix.set_lowerbound(3,0.053843);
  // amatrix.set_uperbound(3,0.253843);
  // amatrix.set_lowerbound(4,-0.20000);
  // amatrix.set_uperbound(4,-0.00000);



  // amatrix.set_lowerbound(0,-0.02000);
  // amatrix.set_uperbound(0,0.00000);
  // amatrix.set_lowerbound(1,0.0);
  // amatrix.set_uperbound(1,0.02);
  // amatrix.set_lowerbound(2,-0.0200);
  // amatrix.set_uperbound(2,-0.000);
  // amatrix.set_lowerbound(3,0.0);
  // amatrix.set_uperbound(3,0.02);
  // amatrix.set_lowerbound(4,-0.220000);
  // amatrix.set_uperbound(4,-0.20000);


  amatrix.set_lowerbound(0,-0.1000);
  amatrix.set_uperbound(0,-0.10000);
  amatrix.set_lowerbound(1,0.127465);
  amatrix.set_uperbound(1,0.127465);
  amatrix.set_lowerbound(2,-0.100703);
  amatrix.set_uperbound(2,-0.100703);
  amatrix.set_lowerbound(3,0.153843);
  amatrix.set_uperbound(3,0.153843);
  amatrix.set_lowerbound(4,-0.10000);
  amatrix.set_uperbound(4,-0.10000);
  
  double **coffs = amatrix.get_coff();
  double **bias = amatrix.get_bias();
  LogicOP *lop =  amatrix.get_lop();
  int num = amatrix.get_num();
  std::cout<<"num of cons is "<<num<<std::endl;
  std::cout<<"before building an absele"<<std::endl;
  AbsEle *aele = new AbsEle(inputDim,name,BOX,num,coffs,lop,bias);
  //  std::cout<<"after building an absele, the ales is "<<*aele<<std::endl<<std::flush;
  
  aeng->iterateWholeNet(aele);
  std::cout<<"computation finished"<<std::endl;
  //  std::cout<<*aele<<std::endl;
  std::cout<<"unsplit is "<<unsplit<<std::endl;
  std::cout<<"split is "<<split<<std::endl;
  return 0;
}

static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << std::endl
              << "Options:\n"
              << "\t-h,--help  \t\tShow this help message\n"
              << "\t-d,--delta  \t\tDelta value, a integer (delta %%) \n"
              << "\t-t,--type   \t\tRobustness type : AI2 / PLANET\n"
//              << "\t-m,--domain \t\tAbstract domain : BOX / TOPO\n"
              << "\t-n,--net    \t\tInputFile, RLV format file\n"
              << "\t-p,--point  \t\tFix input point: a 28*28 MNIST img in integer format\n"
              << "\t-b,--bound  \t\tOutput file, stores the bounds\n"
              << "\t-s,--summary\t\tSummary file: store the statistics\n"
              << std::endl;
}

int main(int argc, char* argv[])
{
    const rlim_t kStackSize = 128L * 1024L * 1024L;
    struct rlimit rl;
    int result = getrlimit(RLIMIT_STACK, &rl);
    if(result == 0){
        if (rl.rlim_cur < kStackSize){
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if(result!=0){
                std::cout<<"stack limit reset faild"<<std::endl;
            }
        }
    }

    if (argc == 1) {
        std::string img("pic1.txt");
        //  std::string net("test.rlv");
        // std::string net("testNetworkB.rlv");
        std::string net("../caffeprototxt/AI2_MNIST_FNN_1/testNetworkB.rlv");

        //  std::string net("MNISTAI2_new.rlv");
        robust_rlv(img,AI, 0.04, net,TOPO,"bounds.txt","summary.txt");
        //  acas();
        return 0;
    }

    double delta = 0.0;
    ROBUST rob_type = PLANET;
    Domain domain_type = BOX;
    std::string net_file_name;
    std::string img_file_name;
    std::string bound_file_name;
    std::string summary_file_name;


    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        } else if ((arg == "-d") || (arg == "--delta")) {
            if (i + 1 < argc) {
                int d = std::atol(argv[++i]);
//                if (delta > 1) {
                delta = d / 1000.0;
//                }
//                cout << delta << endl;
            } else {
                std::cerr << "--delta option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-t") || (arg == "--type")) {
            if (i + 1 < argc) {
                string type_str = string(argv[++i]);
                //cout << type_str << endl;
                if ( type_str == "PLANET") {
                    rob_type = PLANET;
                }
                else if (type_str == "AI2") {
                    rob_type = AI;
                }
                else {
                    std::cerr << "--type option wrong argument." << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "--type option requires one argument." << std::endl;
                return 1;
            }
//        } else if ((arg == "-m") || (arg == "--domain")) {
//            if (i + 1 < argc) {
//                string type_str = string(argv[++i]);
//                //cout << type_str << endl;
//                if (type_str == "BOX") {
//                    domain_type = BOX;
//                }
//                else if (type_str == "TOPO") {
//                    domain_type == TOPO;
//                }
//                else {
//                    std::cerr << "--domain option wrong argument." << std::endl;
//                    return 1;
//                }
//            } else {
//                std::cerr << "--type option requires one argument." << std::endl;
//                return 1;
//            }
        } else if ((arg == "-n") || (arg == "--net")) {
            if (i + 1 < argc) {
                net_file_name = string(argv[++i]);
            } else {
                std::cerr << "--net option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-p") || (arg == "--point")) {
            if (i + 1 < argc) {
                img_file_name = string(argv[++i]);
            } else {
                std::cerr << "--point option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-b") || (arg == "--bound")) {
            if (i + 1 < argc) {
                bound_file_name = string(argv[++i]);
            } else {
                std::cerr << "--bound option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-s") || (arg == "--summary")) {
            if (i + 1 < argc) {
                summary_file_name = string(argv[++i]);
            } else {
                std::cerr << "--summary option requires one argument." << std::endl;
                return 1;
            }
        } else {
            return 1;
        }
    }

#ifdef BDOMAIN
    domain_type = BOX;
#endif

#ifdef TDOMAIN
    domain_type = TOPO;
#endif


    robust_rlv(img_file_name, rob_type, delta, net_file_name, domain_type,bound_file_name,summary_file_name);

}
