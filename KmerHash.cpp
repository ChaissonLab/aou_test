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
        
        
        if (bucket.capacity() == bucket.size()) {
            bucket.reserve(bucket.size() + 10);
        }
        bucket.emplace_back(reminder_, val);
        return &bucket.back().second;
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
    
    totaleles ++;
    if (thekey == 0)  //initiate
    {
        thekey = 1;
    }
    
    else if (thekey == 1) //find collision
    {
        thekey = 2;
        bucketcollision ++;  //record for number of unique collision events
    }
    
    else //N > 1 + 1
    {
        thekey ++ ;  //collision number
    }
    
    
}


void Kmer_hash::initiate()
{
    //flat all collisions
    Hash2.resize(totaleles + 10, 0);
    Hash2_size.resize(totaleles + 10);
    
    uint32_t largestkey = 0;
    uint32_t hash2index = 1;
    for (size_t hash1 = 0; hash1 < large_prime ; ++hash1)
    {
        uint32_t &thekey = Hash1[hash1];
        uint32_t numelement = 0;
        
        if (thekey == 0 )
        {
            continue;
        }
        else
        {
            numelement = thekey;
            thekey = hash2index;
            largestkey = MAX (largestkey, thekey) ;
            Hash2_size[hash2index] = (uint8_t) MIN( MAX_UINT8 , numelement);     //and record the pointer size, condense index by 2
            hash2index += (uint8_t) MIN( MAX_UINT8 , numelement);           //get next local pointer

        }
        
    }
}


void Kmer_hash::addvalue(const ull kmer_int, uint8* counts, Excess_hash &excess_counts)
{
    uint index = findhash(kmer_int);
    
    if (index != 0)
    {
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
    
    uint32_t redirect = thekey;
    for (uint32_t index = redirect;  index < redirect + Hash2_size[redirect]; ++index)
    {
        if (Hash2[index] == 0)
        {
            Hash2[index] = hash2 ;    //unintiated
            return 1;
        }
        else if (Hash2[index] == hash2)
        {
            return 0;
        }
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

uint Kmer_hash::findhash(const ull kmer_int)
{
    uint32_t hash1 = (uint32_t) (kmer_int % large_prime) ;
    uint32_t hash2 = (uint32_t) (kmer_int / large_prime);
    uint32_t thekey = Hash1[hash1];
    uint *posi = NULL;
    if (thekey == 0)
    {
        return 0;
    }
        
    else             //check redirect
    {
        uint32_t redirect = thekey;
        for (uint32_t index = redirect;  index < redirect + Hash2_size[redirect]; ++index)
        {
            if (Hash2[index] == hash2) return index;
        }
    }
    
    return 0;
}


void Kmer_hash::savehash(const std::string& outputfile, ull totalkmers)
{
    std::string tmpfile = outputfile + ".tmp";

    std::ofstream ofs(tmpfile, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: Cannot open temp file: " << tmpfile << std::endl;
        return;
    }

    uint64_t size = Hash2.size();
    uint64_t hash1_size = large_prime;

    // Save external totalkmers
    ofs.write(reinterpret_cast<const char*>(&totalkmers), sizeof(totalkmers));

    // Save internal stats
    ofs.write(reinterpret_cast<const char*>(&totaleles), sizeof(totaleles));
    ofs.write(reinterpret_cast<const char*>(&bucketcollision), sizeof(bucketcollision));

    // Sizes
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    ofs.write(reinterpret_cast<const char*>(&hash1_size), sizeof(hash1_size));

    // Data
    ofs.write(reinterpret_cast<const char*>(Hash1.get()), sizeof(uint) * hash1_size);
    ofs.write(reinterpret_cast<const char*>(Hash2.data()), sizeof(uint) * size);
    ofs.write(reinterpret_cast<const char*>(Hash2_size.data()), sizeof(uint8_t) * size);

    ofs.close();
    std::rename(tmpfile.c_str(), outputfile.c_str());
}

ull Kmer_hash::loadhash(const std::string& inputfile)
{
    std::ifstream ifs(inputfile, std::ios::binary);
    if (!ifs) {
        std::cerr << "Error: Cannot open file for reading: " << inputfile << std::endl;
        return 0;
    }

    ull totalkmers = 0;

    // Load external totalkmers
    ifs.read(reinterpret_cast<char*>(&totalkmers), sizeof(totalkmers));

    // Load internal stats
    ifs.read(reinterpret_cast<char*>(&totaleles), sizeof(totaleles));
    ifs.read(reinterpret_cast<char*>(&bucketcollision), sizeof(bucketcollision));

    uint64_t size = 0;
    uint64_t hash1_size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    ifs.read(reinterpret_cast<char*>(&hash1_size), sizeof(hash1_size));

    if (hash1_size != large_prime) {
        std::cerr << "[loadhash] ERROR: matrix .bin file mismatch." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    Hash1 = std::make_unique<uint[]>(hash1_size);
    Hash2.resize(size);
    Hash2_size.resize(size);

    ifs.read(reinterpret_cast<char*>(Hash1.get()), sizeof(uint) * hash1_size);
    ifs.read(reinterpret_cast<char*>(Hash2.data()), sizeof(uint) * size);
    ifs.read(reinterpret_cast<char*>(Hash2_size.data()), sizeof(uint8_t) * size);

    ifs.close();

    std::cout << "[loadhash] Successfully loaded "
              << "Hash1[" << hash1_size << "], Hash2[" << size << "], "
              << "totalkmers = " << totalkmers << std::endl;

    return totalkmers;
}
