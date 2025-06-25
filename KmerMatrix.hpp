//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#ifndef KmerMatrix_hpp
#define KmerMatrix_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <tuple>
#include <mutex>
#include <memory>
#include <numeric>
using namespace std;


#include "config.hpp"
#include "TreeRound.hpp"

class KmerMatrix
{
    
public:
    
    KmerMatrix()
    {};
    ~KmerMatrix()
    {};
    
    void getNormflat(const uint16* kmervec, const uint16* kmermatrix,const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T  &total_lambda, const node* tree , const uint16 nodenum);
    
    void getNorm(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T &total_lambda);
    
    void WindowCovers(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, const int genenum, const int* results, vector<vector<pair<int,int>>>& covers, const int window);
    
private:
        
    unique_ptr<FLOAT_T[]> row_offsites = unique_ptr<FLOAT_T []>( new FLOAT_T [MAX_UINT16] ) ;
    unique_ptr<FLOAT_T[]> diag_offsites = unique_ptr<FLOAT_T[]>( new FLOAT_T [MAX_UINT16] ) ;
    unique_ptr<uint16[]> repeat_counts = unique_ptr<uint16[]>( new uint16 [MAX_UINT16] ) ;
    unique_ptr<uint16[]> repeat_keys = unique_ptr<uint16[]>( new uint16 [MAX_UINT16] ) ;
    
};

#endif /* KmerMatrix_hpp */
