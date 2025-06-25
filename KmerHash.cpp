//
//  KmerHash.cpp
//  KmerHasher
//
//  Created by walfred on 12/15/24.
//

#include "KmerHash.hpp"

using namespace std;

static uint16_t* Search(std::vector<item> &bucket, const uint8 x)
{
    for (uint i = 0 ; i < bucket.size(); ++i)
    {
        if (bucket[i].first == x ) return &bucket[i].second;
    }
    return NULL;
}

uint16* Excess_hash::add(const uint kmer_int, const uint16_t val )
{
    uint32_t hash_ = kmer_int % modsize;
    uint8 reminder_ = uint8(kmer_int / modsize);

    auto &bucket = buckets[hash_];
    uint16_t* loc = Search(bucket, reminder_);

    if ( loc == NULL)
    {
        totalsize ++;
        
        auto size = bucket.size();
        if (size % 10 == 0) bucket.resize(bucket.size() + 10);
        
        bucket[size] = make_pair(reminder_, val);
                    
        return &bucket[size].second;
    }
    else
    {
        *loc = *loc+val;
        return loc;
    }
}

uint16 Excess_hash::findvalue(const uint kmer_int)
{
    uint16* loc = find(kmer_int);
    
    if (loc == NULL) return 0;
    
    return *loc;
}


uint16* Excess_hash::find(const uint kmer_int)
{
    uint hash_ = kmer_int % modsize;
    uint8 reminder_ = uint8(kmer_int / modsize);
    
    auto &bucket = buckets[hash_];
    uint16* loc = Search(bucket, reminder_);
    
    return loc;
}



void Kmer_hash::preview(const ull kmer_int)
{
    uint32_t hash1 = (uint32_t) (kmer_int % large_prime) ;
    uint32_t hash2 = (uint32_t) (kmer_int / large_prime);
    uint32_t &thekey = Hash1[hash1];
    

    if (thekey == MAX_UINT32)  //initiate
    {
        thekey = hash2;
    }
    
    else if (thekey < 0x80000000) //find collision
    {
        thekey = 0x80000002;
        bucketcollision ++;  //record for number of unique collision events
        totalcollision += 2;  //record for number of total collision
    }
    
    else //N > 1 + 1
    {
        thekey ++ ;  //collision number
        totalcollision ++;
    }
    
    
}


void Kmer_hash::initiate()
{
    //flat all collisions
    Hash2.resize(totalcollision + 10, MAX_UINT32);
    Hash2_size.resize(totalcollision/2 + 10);
    
    uint32_t largestkey = 0;
    uint32_t hash2index = 0;
    for (size_t hash1 = 0; hash1 < large_prime ; ++hash1)
    {
        uint32_t &thekey = Hash1[hash1];
        uint32_t numelement = 0;
        
        if (thekey == MAX_UINT32 )
        {
            continue;
        }
        else if (thekey > 0x80000000)
        {
            numelement = thekey - 0x80000000;
            thekey = hash2index + 0x80000000;
            largestkey = MAX (largestkey, thekey) ;
            Hash2_size[hash2index/2] = (uint8_t) MIN( MAX_UINT8 , numelement);     //and record the pointer size, condense index by 2
            hash2index += (uint8_t) MIN( MAX_UINT8 , numelement);           //get next local pointer

        }
    }
}

uint Kmer_hash::findhash(ull kmer_int)
{
    uint * loc = find(kmer_int);
    
    if (loc == NULL) return 0;
    
    return *loc;
}

void Kmer_hash::addvalue(const ull kmer_int, uint8* counts, Excess_hash &excess_counts)
{
    uint index = findhash(kmer_int);
    
    if (index != NULL)
    {
        if ( index == 0) return;
        
        if ( __builtin_expect(counts[index] < 255, 1))
        {
            counts[index] ++;
        }
        else
        {
            excess_counts.add(index, 1);
        }
    }
}

bool Kmer_hash::add(const ull kmer_int)
{
    uint32_t hash1 = (uint32_t) (kmer_int % large_prime) ;
    uint32_t hash2 = (uint32_t) (kmer_int / large_prime);
    uint32_t& thekey = Hash1[hash1];
    
    if (thekey != MAX_UINT32)      //redirected
    {
        if (thekey < 0x80000000)
        {
            assert(thekey == hash2);
            return 0;
        }
        
        uint32_t redirect = thekey - 0x80000000;
        for (uint32_t index = redirect;  index < redirect + Hash2_size[redirect/2]; ++index)
        {
            if (Hash2[index] == MAX_UINT32)
            {
                Hash2[index] = hash2 ;    //unintiated
                return 1;
            }
            else if (Hash2[index] == thekey)
            {
                return 0;
            }
        }
        
    }
    else
    {
        thekey = hash2;
        return 1;
    }
    
    return 0;
}


uint16 Kmer_hash::findvalue(const ull kmer_int, const uint8* counts, Excess_hash &excess_counts)
{
    uint loc = findhash(kmer_int);
    
    uint8 count = counts[loc];
    
    if ( __builtin_expect(count < 255, 1))
    {
        return count;
    }
    else
    {
        return count + excess_counts.findvalue(loc);
    }
    
    return count;
}

uint *Kmer_hash::find(const ull kmer_int)
{
    uint32_t hash1 = (uint32_t) (kmer_int % large_prime) ;
    uint32_t hash2 = (uint32_t) (kmer_int / large_prime);
    uint32_t thekey = Hash1[hash1];
    uint *posi = NULL;
    if (thekey == MAX_UINT32)
    {
        return NULL;
    }
        
    else if (thekey >= 0x80000000)             //check redirect
    {

        uint32_t redirect = thekey - 0x80000000;
        for (uint32_t index = redirect;  index < redirect + Hash2_size[redirect/2]; ++index)
        {
            if (Hash2[index] == hash2) return &Hash2[index];
        }
    }
    
    else
    {
        if (hash2 == Hash1[hash1])
        {
            return &Hash1[hash1];
        }

    }
    
    return posi;
}
