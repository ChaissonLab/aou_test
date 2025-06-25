//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#ifndef Regression_hpp
#define Regression_hpp

#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "config.hpp"

using namespace std;

class Regression
{
public:
    Regression()
    {};
    ~Regression()
    {};
    
    int lawson_hanson_nnls(const FLOAT_T *kernal_vec, const FLOAT_T *weightnorm, uint16 size, FLOAT_T *coefs, FLOAT_T *residuel);

    void ruduce_nnls(FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T* reduce_matrix, const uint16 gnum, FLOAT_T * coefs, FLOAT_T * residuels, const uint16 numsmallgroups, const vector<uint16> &smallgroups);
    
    void Call(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T* reduce_matrix, FLOAT_T &total_lambda, const uint *kmercounts, FLOAT_T * coefs, FLOAT_T * residuels, const size_t numgroup, const vector<uint16> &groups, const uint16 numsmallgroups, const vector<uint16> &smallgroups, const vector<uint> &groupkmernums);

private:
    
    vector<bool> active_or_passive = vector<bool> (MAX_UINT16);
        
    vector<uint16> passive_set = vector<uint16> (MAX_UINT16);
    
    vector<uint16> passive_set_old = vector<uint16> (MAX_UINT16);
    
    vector<uint16> smallgroup_count = vector<uint16> (MAX_UINT16);
    
    std::unique_ptr<FLOAT_T[]> reduce_vec = std::make_unique<FLOAT_T[]>(MAX_UINT16);

    std::unique_ptr<FLOAT_T[]> reduce_coef = std::make_unique<FLOAT_T[]>(MAX_UINT16);
};



#endif /* Regression_hpp */
