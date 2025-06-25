//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#ifndef KmerCounter_hpp
#define KmerCounter_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <tuple>
#include <filesystem>


#if defined(__APPLE__) && defined(__clang__)
// Apple Clang uses __fs::filesystem instead of std::filesystem
namespace std {
    namespace filesystem = __fs::filesystem;
}
#endif


#include "FastaReader.hpp"
#include "FastqReader.hpp"
#include "CramReader.hpp"
#include "KtableReader.hpp"
#include "KmerHash.hpp"

struct Hash10M {
    std::size_t operator()(ull key) const {
        return static_cast<std::size_t>(key % prime10M);
    }
};

template <typename T1>
string int_to_kmer(T1 value, int size = 31)
{
    string kmer (31,' ');
    for (int i = 0; i<31 ; ++i)
    {
        kmer[30-i] = "ACGT"[value % 4];
        value /= 4;
    }
    return kmer;
}

static inline bool base_to_int(const char base, int &converted)
{
    
    switch (base)
    {
        case 'A' :
            converted=0b00;
            break;
        case 'C' :
            converted=0b01;
            break;
        case 'G' :
            converted=0b10;
            break;
        case 'T' :
            converted=0b11;
            break;
        default:
            return 0;
    }
    
    return 1;
}

class KmerCounter
{
public:
    
    KmerCounter()
    {};
    ~KmerCounter()
    {};
    
    ull read_target(KtableReader &fastafile);
    
    ull read_target(FastaReader &fastafile);
    
    ull read_target(const char* infile);
	
    void load_backgrounds(const char * backfile);
        
    void count_kmer(CramReader &file, counterint* samplevecs, Excess_hash &exbucks , ull_atom &nBases, ull_atom &nReads ,ull_atom &nBg);
    
    void count_kmer(FastaReader &file, counterint* samplevecs, Excess_hash &exbucks, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg);
   
    void count_kmer(FastqReader &file, counterint* samplevecs, Excess_hash &exbucks , ull_atom &nBases, ull_atom &nReads ,ull_atom &nBg);
 
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes,std::vector<float>& deps,int numthread);
    
    void Call(const char* infile, counterint* samplevecs, Excess_hash &exbucks, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads);
        
    std::vector<ull> backgrounds;
	ull totalbks = 0;
    uint totalkmers = 0;
    std::vector<char *> regions;
    
    const int klen = 31;
    
    Kmer_hash kmer_hash;
    std::mutex counting_lock;
    
    template <class typefile>
    void count_kmer_(typefile &file, counterint* samplevecs, Excess_hash &exbucks, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg,  int nthreads);

};


    





#endif /* KmerCounter_hpp */
