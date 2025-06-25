//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#include <Eigen/Dense>
#include <vector>
#include <limits>
#include <algorithm>
#include "Regression.hpp"
#include <iomanip>
#include <functional>
#include <chrono>
#include <thread>

#define max_repetitions 5
extern bool optioncorr;


template<typename Func>
void try_function_with_retries(Func func, int max_attempts = 100)
{
    int attempts = max_attempts;
    
    while (attempts-- > 0)
    {
        try
        {
            func();
            return;
        }
        
        catch (const std::bad_alloc&)
        {
            std::cerr << "Matrix solving memory allocation failed in function call, Attempt: " << max_attempts - attempts << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(10));
        }
    }
    
    std::cerr << "Matrix solving failed to call function after " << max_attempts << " attempts. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
}


void trial_solution(vector<uint16>& passive_set, const uint16 passive_num, const Vector_T &y, const Matrix_T &A, const uint16 size, Vector_T &s)
{
    assert( passive_num <= size );
    
    s.head(size).setZero();
    
    
    if (passive_num == 1)
    {
        s(passive_set[0]) = (FLOAT_T)y(passive_set[0]) / (FLOAT_T)A(passive_set[0], passive_set[0]);
        return ;
    }
     
    
    Matrix_T sub_matrix = Matrix_T(passive_num, passive_num);
    Vector_T sub_vector = Vector_T(passive_num);
    
    for (uint16 i = 0; i < passive_num; ++i)
    {
        sub_vector(i) = y(passive_set[i]);
        
        for (uint16 j = 0; j < passive_num; ++j)
        {
            sub_matrix(i, j) = A(passive_set[i], passive_set[j]);
        }
    }
    
    Vector_T s_ = sub_matrix.colPivHouseholderQr().solve(sub_vector);
        
    for (size_t i = 0; i < passive_num; ++i)
    {
        s(passive_set[i]) = s_(i);
    }
}

int Regression::lawson_hanson_nnls(const FLOAT_T *kernal_vec, const FLOAT_T *weightnorm, uint16 size, FLOAT_T *coefs, FLOAT_T *residuel)
{
    int max_iterations = MIN(2*size, 1000);
    
    const FLOAT_T tol = size*numeric_limits<FLOAT_T>::epsilon();
    
    const Vector_T y = Eigen::Map<Vector_T> ((FLOAT_T*) kernal_vec, size);
    
    const Matrix_T A = Eigen::Map<Matrix_T> ((FLOAT_T*) weightnorm, size, size);
        
    Vector_T x = Vector_T::Zero(size), x_trial(size) , r (size);
    
    uint16 no_update = 0;
    
    active_or_passive.assign(size,0); //0 for active, 1 for passive
    uint16 passive_num = 0;
    uint16 passive_num_old = 0;
    
    r = y - A * x;
    
    int j = 0;
    while (passive_num < size && j++ < max_iterations)
    {
        // Find the max residual max_r in active set and its index;
        
        int max_r_index = 0;
        FLOAT_T max_r = std::numeric_limits<FLOAT_T>::min();

        for (uint16 i = 0; i < size; ++i)
        {
            if (active_or_passive[i]) continue; //in passive, ignore
            
            if (r[i] > max_r)
            {
                max_r = r[i];
                max_r_index = i;
            }
        }
                
        if (max_r <= tol)  break; //non-negative achieved

        //make a copy for passive_set for late comparison
        passive_num_old = passive_num;
        passive_set_old.assign(passive_set.begin(), passive_set.begin() + passive_num_old);
        
        // Move the variable j from the passive set to active set .
        passive_set[passive_num++] = max_r_index;
        active_or_passive[max_r_index] = 1;
        
        
        try_function_with_retries( [&]() {trial_solution(passive_set, passive_num, y, A, size, x_trial); } );
        
        //trial_solution(passive_set, passive_num, y, A, size, x_trial); //solve || y - A * x || on passive vectors
        
        int k = 0;
        while (passive_num && k++ < max_iterations)
        {

            // Find the index i in the passive set P that has the smallest alpha value.
            FLOAT_T alpha = std::numeric_limits<FLOAT_T>::max();
            int alpha_index = -1;
            for (size_t i = 0; i < passive_num; ++i)
            {
                int idx = passive_set[i];
                FLOAT_T alpha_i = x_trial[idx] ;
                if (alpha_i < alpha)
                {
                    alpha = alpha_i;
                    alpha_index = idx;
                }
            }
            
                        
            if (alpha > tol)
            {
                break;
            }
                        
            int passive_shrink_index = 0;
            for (size_t i = 0; i < passive_num; ++i)
            {
                int idx = passive_set[i];
                if (x_trial[idx] > tol)
                {
                    passive_set[passive_shrink_index++] = passive_set[i];
                }
                else
                {
                    active_or_passive[max_r_index] = 0;
                }
            }
            
            passive_num = passive_shrink_index;
            
            // Update x.
            x = x + alpha * (x_trial - x);
            
            if (passive_num) try_function_with_retries( [&]() {trial_solution(passive_set, passive_num, y, A, size, x_trial); } );
            
        }

        x = x_trial.eval();
        // Update the residual.
        r = y - A * x;
        
        if(passive_num == passive_num_old  && equal(passive_set_old.begin(), passive_set_old.begin() + passive_num, passive_set.begin() ) )
        {
            no_update ++;
        }
        else
        {
            no_update = 0;
        }
        
        
        if (no_update > max_repetitions) break;
        
    }
    
    r = y - A * x;
    
    for (int i = 0; i < size; ++i)
    {
        coefs[i] = x(i);
        residuel[i] = r(i);
    }
    
    
    return 0;
}


inline void getlocalnum_mul(const FLOAT_T* coefs, const uint16* rowdata,const FLOAT_T totalnum, FLOAT_T &localnum, const vector<FLOAT_T>& grouptotalnums, vector<FLOAT_T>& grouplocalnums, const vector<uint16> &groups)
{
    
    uint16 geneindex = 0;
    
    switch (rowdata[0])
    {
        case '_': case '=':
            break;

        case '-':
                        
            localnum = totalnum;
            grouplocalnums = grouptotalnums;
            
            for (int i = rowdata[9] ; i < rowdata[1] ; i++)
            {
                geneindex = rowdata[ FIXCOL + i];
                localnum -= coefs[geneindex];
                grouplocalnums[groups[geneindex]] -= coefs[geneindex];
            }
            
            break;
            
        case '+':
            
            localnum = 0;
            std::fill(grouplocalnums.begin(), grouplocalnums.end(), 0);

            for (int i = rowdata[9]  ; i < rowdata[1] ; i++)
            {
                geneindex = rowdata[ FIXCOL + i];
                localnum += coefs[geneindex];
                grouplocalnums[groups[geneindex]] += coefs[geneindex];
            }
            
            break;
        default:
            
            break;
    }
}

inline void GetMedianNofil_mul(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, const FLOAT_T &totalnum, const vector<FLOAT_T> &grouptotalnums,const vector<uint16> &groups, const uint16 groupnum, vector<size_t> &grouptotalobs, vector<vector<FLOAT_T>> &allratios)
{
    
    FLOAT_T localnum = 0.0;
    vector<FLOAT_T> grouplocalnums (groupnum + 1, 0.0);
    
    for (size_t i = 0; i < knum; ++i)
    {
        getlocalnum_mul(coefs, rowdata, totalnum, localnum, grouptotalnums, grouplocalnums, groups);
        
        if (localnum >= 0.5 && kmervec[i] >= 2)
        {
            grouplocalnums[groupnum] = localnum;
            
            FLOAT_T ratio = kmervec[i] / localnum;
            
            if (optioncorr && (rowdata[5] & 0x3F) >= errorcutoff1 )
            {
                float corr = 0.01 * (rowdata[5] & 0xFFC0)/64;
                ratio *= corr;
            }

            for (int i = 0 ; i < groupnum + 1; ++i)
            {
                if (grouplocalnums[i] >= 0.5 && grouptotalnums[i] >= corrstartpoint)
                {
                    allratios[i][grouptotalobs[i]++] = ratio;
                }
            }
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
    
}

void resize_mul(vector<vector<FLOAT_T>> &allratios, const vector<uint> &groupkmermax, const vector<FLOAT_T> &grouptotalnums, const uint knum, const uint16 groupnum)
{
    for (uint16 index = 0; index < groupnum;++index)
    {
        if ( grouptotalnums[index] >= corrstartpoint) allratios[index].resize(10 * sufficient + groupkmermax[index],0);
        assert(groupkmermax[index] <= knum);
    }
    allratios[groupnum].resize(10 * sufficient + knum,0);
}

void aggregateCorr_mul(FLOAT_T * coefs, const uint16* kmervec, const uint16* kmermatrix, const vector<uint> &groupkmernums, const uint16 gnum, const uint knum, const vector<uint16> &groups, const uint16 groupnum, vector<FLOAT_T> &corrections)
{
    
    vector<size_t> grouptotalobs (groupnum + 1, 0);
    vector<FLOAT_T> grouptotalnums (groupnum + 1, 0.0);
    vector<vector<FLOAT_T>> allratios(groupnum + 1);
    
    FLOAT_T totalnum = 0.0;
    for (int i = 0; i < gnum; ++i)
    {
        totalnum += coefs[i];
        grouptotalnums[groups[i]] += coefs[i];
    }
    grouptotalnums[groupnum] = MAX(totalnum, 1.0);
    
    //try_function_with_retries( [&]() { resize_mul(allratios, groupkmernums, grouptotalnums, knum,groupnum);} );

    resize_mul(allratios, groupkmernums, grouptotalnums, knum,groupnum);
    
    GetMedianNofil_mul(coefs, kmervec,kmermatrix,knum, totalnum,grouptotalnums,groups, groupnum, grouptotalobs, allratios);
    
    for (uint16 i = 0; i < groupnum + 1; ++i)
    {
        if (grouptotalnums[i] >= corrstartpoint)
        {
            std::sort(allratios[i].begin(), allratios[i].begin() + grouptotalobs[i]);
            
            corrections[i] = allratios[i][grouptotalobs[i]/2];
        }
    }
}

inline void getlocalnum(const FLOAT_T* coefs, const uint16* rowdata,const FLOAT_T totalnum, FLOAT_T &localnum)
{
    switch (rowdata[0])
    {
        case '_': case '=':
            break;

        case '-':
            
            localnum = totalnum;
            for (int i = rowdata[9] ; i < rowdata[1] ; i++)
            {
                localnum -= coefs[rowdata[ FIXCOL + i]];
            }
            
            break;
            
        case '+':
            
            localnum = 0;
            for (int i = rowdata[9] ; i < rowdata[1] ; i++)
            {
                localnum += coefs[rowdata[ FIXCOL + i]];
            }
            
            break;
        default:
            
            break;
    }
}


inline void GetMedianNofil(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, FLOAT_T &totalnum,size_t &totalobs, vector<FLOAT_T>& ratios)
{
    FLOAT_T localnum = 0.0;
    for (size_t i = 0; i < knum; ++i)
    {
        
        getlocalnum(coefs, rowdata, totalnum, localnum);
        
        if ( localnum >= 0.5 && kmervec[i] >= 2)
        {
            FLOAT_T ratio = kmervec[i] / localnum;
            
            if (optioncorr && (rowdata[5] & 0x3F) >= errorcutoff1 )
            {
                float corr = 0.01 * (rowdata[5] & 0xFFC0)/64;
                ratio *= corr;
            }
                
            ratios[totalobs] = ratio;
            totalobs ++;
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
    
    
}

FLOAT_T aggregateCorr(const FLOAT_T * coefs, const uint16* kmervec, const uint16* kmermatrix, const uint16 gnum, const uint knum)
{
    
    vector<FLOAT_T> ratios(knum, 0.0);
        
    size_t totalobs  = 0;
    FLOAT_T totalnum = 0.0;
    for (int i = 0; i < gnum; ++i)
    {
        totalnum += coefs[i];
    }
    
    GetMedianNofil(coefs, kmervec, kmermatrix, knum, totalnum, totalobs, ratios);
    
    std::sort(ratios.begin(), ratios.begin() + totalobs);
    
    
    FLOAT_T median = ratios[totalobs/2];
    

    return median;
}


void Regression::ruduce_nnls( FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T* reduce_matrix, const uint16 gnum, FLOAT_T * coefs, FLOAT_T * residuels, const uint16 numsmallgroups, const vector<uint16> &smallgroups)
{
    
    std::fill(smallgroup_count.begin(), smallgroup_count.begin() + numsmallgroups, 0);
    std::fill_n(reduce_vec.get(), numsmallgroups, static_cast<FLOAT_T>(0));
    std::fill_n(reduce_coef.get(), numsmallgroups, static_cast<FLOAT_T>(0));
    
    for (auto index: smallgroups)
    {
        smallgroup_count[index] ++;
    }
    
    for (int i = 0; i < gnum ; ++i)
    {
        auto i_smallgroup = smallgroups[i];
        reduce_vec[i_smallgroup] += norm_vec[i];
        for (int j = i; j < gnum ; ++j)
        {
            auto j_smallgroup = smallgroups[j];
            
            reduce_matrix[i_smallgroup * numsmallgroups + j_smallgroup] += norm_matrix[i * gnum + j];
            
        }
    }
    
    for (int i = 0; i < numsmallgroups ; ++i)
    {
        reduce_vec[i] /= smallgroup_count[i];
        
        for (int j = i; j < numsmallgroups ; ++j)
        {
            reduce_matrix[i * numsmallgroups + j] /= (smallgroup_count[i] * smallgroup_count[j]);
            reduce_matrix[j * numsmallgroups + i] = reduce_matrix[i * numsmallgroups + j];
        }
    }
    
    
    lawson_hanson_nnls(reduce_vec.get(), reduce_matrix, numsmallgroups, reduce_coef.get(), residuels);
    
    std::fill_n(residuels, gnum, static_cast<FLOAT_T>(0));
    
    for (int i = 0; i < gnum ; ++i)
    {
        auto i_smallgroup = smallgroups[i];
        coefs[i] = reduce_coef.get()[i_smallgroup] / smallgroup_count[i_smallgroup];
        
        for (int j = 0; j < gnum ; ++j)
        {
            residuels[j] += coefs[i] * norm_matrix[gnum * i + j];
        }
    }
    for (int i = 0; i < gnum ; ++i)
    {
        residuels[i] -=  norm_vec[i];
    }
}

void Regression::Call(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T* reduce_matrix,FLOAT_T &total_lambda, const uint *kmercounts, FLOAT_T * coefs, FLOAT_T * residuels, const size_t numgroup, const vector<uint16> &groups, const uint16 numsmallgroups, const vector<uint16> &smallgroups,const vector<uint> &groupkmernums)
{
    
    float regressed_kmer = 0.0;
    
    if (gnum > 10000)
    {
        ruduce_nnls(norm_vec, norm_matrix, reduce_matrix, gnum, coefs, residuels, numsmallgroups, smallgroups);
    }
    else
    {
        lawson_hanson_nnls(norm_vec, norm_matrix, gnum, coefs, residuels);
    }
    
    auto t1 = std::chrono::steady_clock::now();
    /*
    float regressed_kmer = 0.0;
    
    for (int i = 0; i < gnum; ++i)
    {
        regressed_kmer += kmercounts[i] * coefs[i];
    }
    regressed_kmer = (regressed_kmer > 0 ) ? regressed_kmer : 1;
     
    float correction = total_lambda/regressed_kmer;
    */
    if (numgroup <= 1)
    {
        FLOAT_T correction = ( aggregateCorr(coefs, kmervec, kmermatrix, gnum, knum) ) / depth;
        
        for (int i = 0; i < gnum; ++i)
        {
            
            coefs[i] *= correction;
            
            residuels[i] *= correction;
            
        }
    }
    
    
    else
    {
        vector<FLOAT_T> corrections(numgroup + 1, 0.0);
        
        aggregateCorr_mul(coefs, kmervec, kmermatrix, groupkmernums, gnum, knum, groups, numgroup, corrections);
        
        for (int index = 0; index < numgroup ; ++index)
        {
            corrections[index] = (corrections[index]  )/depth;
        }
        
        FLOAT_T correction = ( corrections[numgroup] )/depth ;
        
        
        for (int i = 0; i < gnum; ++i)
        {
            
            if (corrections[groups[i]] > 0)
            {
                coefs[i] *= corrections[groups[i]];
                
                residuels[i] *= corrections[groups[i]];
            }
            else
            {
                coefs[i] *= correction;
                
                residuels[i] *= correction;
            }
        }
        
    }
    
    
    auto t3 = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t1).count();
    std::cout << "Correctiong for copy number: "<< " at " << elapsed << " ms" << std::endl;
    
}

