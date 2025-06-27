//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License.
//  If you use this code, please cite our work.
//

#ifndef Processor_hpp
#define Processor_hpp

#include <atomic>
#include <mutex>
#include <thread>
#include <atomic>
#include <memory>
#include <stdio.h>
#include <string>
#include <vector>
#include <chrono>

#include "config.hpp"
#include "FastaReader.hpp"
#include "FastqReader.hpp"
#include "KtableReader.hpp"
#include "KmerMatrix.hpp"
#include "KmerCounter.hpp"
#include "KmerWindow.hpp"
#include "Regression.hpp"
#include "PriorData.hpp"
#include "Genotyper.hpp"

using namespace std;


class Processor
{
public:
    Processor
    (
              std::vector<std::string>& infiles,
              std::vector<std::string>& outfiles,
              std::vector<float> &d,
              std::string &mfile,
              std::string &bfile,
              std::unordered_set<std::string> &g,
              std::vector<char *> &r,
              const int w,
              const int n,
              const int N,
	      string ref=""
    ):
    inputfiles(infiles),
    outputfiles(outfiles),
    depths(d),
    genes(g),
    matrixfile(mfile),
    backgroundfile(bfile),
    priordata_manager(mfile, (n+1)/2),
    regions(r),
    window(w),
    nthreads(n),
    Nsubthreads(N),
    reference(ref)
    {};
    
    
    void Run();
    void Load();
    void Onethread();
    
    const std::vector<std::string>& inputfiles;
    const std::vector<std::string>& outputfiles;
    const std::string &matrixfile;
    const std::string &backgroundfile;
    std::vector<char *> &regions;
    const std::vector<float> &depths;
    const std::unordered_set<std::string> &genes;
    const int window;
    const int nthreads;
    const int Nsubthreads;
  std::string reference;
private:
    uint totalkmers, totalgroups;
    std::atomic_uint restfileindex = {0};
    std::mutex Threads_lock;
    std::mutex Threads_lock2;
    
    //kmer_hash_type kmer_hash;
    //kmer_hash_type_mul kmer_multi_hash;
    //unordered_set<ull, Hash10M> backgroud;
    KmerCounter *Counter =NULL;
    PriorData priordata_manager;
};





#endif /* Processor_hpp */

