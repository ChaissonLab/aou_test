//
//  KmerHash.hpp
//  KmerHasher
//
//  Created by walfred on 12/15/24.
//

#ifndef KmerHash_hpp
#define KmerHash_hpp

#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <functional>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <iostream>


#include "config.hpp"

#define MAX_UINT22 4194303
#define uint unsigned int

using namespace std;

struct uint40
{
    uint40() : first(0), second(0) {}
    uint40 (const ull x): first( x % UINT_MAX ) , second( x / UINT_MAX )
    {};
    uint40 (const uint32_t x, const uint8_t y): first( x ) , second(y)
    {};
    
    auto operator==(const uint40& other) const
    {
        return (this->first == other.first && this->second == other.second);
    }
 
    uint first;
    uint8 second;
};

typedef std::pair<uint8,uint16> item;

class Excess_hash
{
public:
    Excess_hash() : modsize(bucketprime)
    {
        buckets.resize(modsize);
    };

    ~Excess_hash() = default; // unique_ptr handles cleanup

    uint16* add(const uint kmer_int, const uint16_t index);
    uint16* find(const uint kmer_int);
    uint16 findvalue(const uint kmer_int);
    void reset(){buckets.clear();buckets.resize(modsize);};

    const size_t modsize;

    std::vector<std::vector<item>> buckets;

    uint32_t totalsize = 0;
};

class Kmer_hash
{
public:
    Kmer_hash()
        : Hash1(make_unique<uint[]>(large_prime))
    {
        std::fill_n(Hash1.get(), large_prime, 0);
    };
    ~Kmer_hash()
    {
        
    };

	void savehash(const std::string& outputfile, ull totalkmers);
    ull loadhash(const std::string& outputfile);
	void addvalue(const ull kmer_int, uint8* counts, Excess_hash &excess_counts);
    uint16 findvalue(const ull kmer_int, const uint8* counts, Excess_hash &excess_counts);
    void preview(const ull kmer_int);
    void initiate();
    bool add(const ull kmer_int);
	uint findhash(ull kmer_int);
    ull totalkmer = 0;
    uint32_t totaleles = 0;
    uint32_t bucketcollision = 0;
    unique_ptr<uint[]> Hash1;
    vector<uint> Hash2;
    vector<uint8_t> Hash2_size;
    
};




#endif /* KmerHash_hpp */

