//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#include "KmerMatrix.hpp"

#include <fstream>
#include <string>
#include <unordered_set>

extern bool optioncorr;

inline void getEachRowValue(const FLOAT_T depth, const int count, const char sign, const uint16 uniqcounter,const uint16 flag,const FLOAT_T mean_repeat,FLOAT_T &total_lambda, FLOAT_T &norm_value, FLOAT_T &weight_value)
{
    float ori_weight = 1.0;
    if (sign==1 && uniqcounter == 1)
    {
        ori_weight = 0.05 ;
    }
    
    if (count < 3 )
    {
        if (mean_repeat != 2.0) weight_value += (4 * ori_weight/mean_repeat - ori_weight);
        return ;
    }
    
    
    float count_f;           //float copy number value
    float count_i;   //estimate number of copy and at least one copy
    float new_weight;
    
    count_f = 1.0 * count/depth  ;
    
    if (optioncorr && (flag & 0x3F) > errorcutoff1 )
    {
        
        float corr = 0.01 * (flag & 0xFFC0)/64;
        count_f *= corr;
    }
    
    count_i = (int(count_f+ 0.5) > 1.0) ? int(count_f+ 0.5) : 1.0;   //estimate number of copy and at least one copy
    
    new_weight =  1.00/(count_i*count_i * mean_repeat);
    
    norm_value += 4 * count_i * new_weight ;
    weight_value += 4 *  new_weight - ori_weight;
    //weight_value += (new_weight - ori_weight);  //weight is reversely proportion to square of estimated copy number, we calculate offsite to the original weight
               //norm vector value of this kmer = count_i * weight = 1.00/count_i
        
    total_lambda += count_f;
}

//A function to change values of norm vec and norm matrix for a list of kmers found in exactly the same list of samples
//This version is for binaray major kmers (frequencies > 50% and no sample has more than 1)
//We reverse calculate values for samples that missing this kmer
inline void getEachRowNorm_major(const uint16 rowsize, const uint16 *row, const uint16 gnum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, const FLOAT_T norm_value, const FLOAT_T  weight_value, FLOAT_T  &vec_offsite, FLOAT_T&matrix_offsite, FLOAT_T* row_offsites)
{
    
    vec_offsite += norm_value;            //major kmer, default is add this weight change to whole norm vector
    for (int i = 0; i < rowsize; ++i)
    {
        norm_vec[row[i]] -= norm_value;
    }
    
    if (weight_value == 0)
    {
        return;
    }
    
    matrix_offsite += weight_value;            //major kmer, default is add this weight change to whole matrix
    
    for (int i = 0; i < rowsize; ++i)            //sample missing, ignored its row/col from the weight change
    {
        norm_matrix[row[i]*gnum+row[i]] -= weight_value * 0.5;   //we will double diagonal elements after
        row_offsites[row[i]] -= weight_value;               //row[i] sample is missing this kmer, should not be affected, use offsite to correct back
        
        for (int j = i+1; j < rowsize; ++j)
        {
            norm_matrix[row[i]*gnum+row[j]] += weight_value;           //we use add to replace matrix multiplication because most kmers only found once in a locus, when row[i] == row[j], it should be doubled, so we will double diagonal elements after
        }
    }
    
}


//This version is for minor or non-binaray major kmers
inline void getEachRowNorm_minor_repeat(const uint16 rowsize, const uint16 *row, const uint16 gnum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, const FLOAT_T norm_value, const FLOAT_T  weight_value, uint16 *repeat_keys, uint16 *repeat_counts)
{
    
    uint unique_count = 0;
    uint16 lastrow = MAX_UINT16;
    for (int i = 0; i < rowsize; ++i)
    {
        if (row[i] != lastrow)
        {
            repeat_counts[unique_count] = 1;
            repeat_keys[unique_count++] = row[i];
        }
        else
        {
            repeat_counts[unique_count] ++;
        }
    }
    
    for (int i = 0; i < unique_count; ++i)
    {
        auto key_i = repeat_keys[i];
        auto count_i = repeat_counts[i];
        norm_matrix[key_i*gnum+key_i] += count_i * count_i * weight_value/2;
        
        for (int j = i+1; j < rowsize; ++j)
        {
            auto key_j = repeat_keys[j];
            auto count_j = repeat_counts[j];
            norm_matrix[key_i*gnum+key_j] += count_i * count_j *weight_value;
        }
    }
    
    
}

//This version is for minor or non-binaray major kmers
inline void getEachRowNorm_minor(const uint16 rowsize, const uint16 *row, const uint16 gnum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, const FLOAT_T norm_value, const FLOAT_T  weight_value, const FLOAT_T mean_repeat, uint16 *repeat_keys, uint16 *repeat_counts)
{
    
    for (int i = 0; i < rowsize; ++i)
    {
        norm_vec[row[i]] += norm_value;
    }
    
    if (weight_value == 0)
    {
        return;
    }
    if (mean_repeat > 3.0)
    {
        getEachRowNorm_minor_repeat(rowsize, row,  gnum, norm_vec,  norm_matrix,norm_value, weight_value, repeat_keys, repeat_counts);
    }
    
    else
    {
        for (int i = 0; i < rowsize; ++i)
        {
            norm_matrix[row[i]*gnum+row[i]] += weight_value/2;
            
            for (int j = i+1; j < rowsize; ++j)
            {
                norm_matrix[row[i]*gnum+row[j]] += weight_value;
            }
        }
    }
}


void getPriorNorm2(const uint numposi, const size_t colsize, FLOAT_T* prior_norm, const uint16* temprow, uint16 matrixsize, FLOAT_T weight)
{
    for (auto index = 0; index < numposi; ++index)
    {
        
        uint16 i = temprow[index];
        prior_norm[matrixsize*i + i] += 0.5*weight ;
        
        for (size_t jndex =index+1; jndex< numposi ; ++jndex)
        {
            uint16 j = temprow[jndex];
            prior_norm[matrixsize*i + j] += weight;
            
        }
    }
    
    for (auto index = numposi; index < colsize; ++index)
    {
        uint16 i = temprow[index];
        prior_norm[matrixsize*i + i] += 0.5*weight ;
        
        for (size_t jndex =index+1; jndex< colsize ; ++jndex)
        {
            uint16 j = temprow[jndex];
            prior_norm[matrixsize*i + j] += weight;
        }
    }
    for (auto index = 0; index < numposi; ++index)
    {
        uint16 i = temprow[index];
        
        for (size_t jndex =numposi; jndex< colsize ; ++jndex)
        {
            uint16 j = temprow[jndex];
            prior_norm[matrixsize* (MIN(i,j)) + (MAX(i,j)) ] -= weight;
            
        }
    }
    
    
}

//This version is for minor or non-binaray major kmers
inline void getEachRowNorm(const uint16 posisize, const uint16 rowsize, const uint16 *row, const uint16 gnum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, const FLOAT_T norm_value, const FLOAT_T  weight_value)
{
    
    for (int i = 0; i < posisize; ++i)
    {
        norm_vec[row[i]] += norm_value;
    }
    
    for (int i = posisize; i < rowsize; ++i)
    {
        norm_vec[row[i]] -= norm_value;
    }
    
    if (weight_value == 0)
    {
        return;
    }
    

    for (int i = 0; i < posisize; ++i)
    {
        norm_matrix[row[i]*gnum+row[i]] += weight_value/2;
        for (int j = i+1; j < posisize; ++j)
        {
            norm_matrix[row[i]*gnum+row[j]] += weight_value;
        }
    }
    
    for (int i = posisize; i < rowsize; ++i)
    {
        norm_matrix[row[i]*gnum+row[i]] += weight_value/2;
        for (int j = i+1; j < rowsize; ++j)
        {
            norm_matrix[row[i]*gnum+row[j]] += weight_value;
        }
    }
    
    for (int i = 0; i < posisize; ++i)
    {
        for (int j = posisize; j < rowsize; ++j)
        {
            norm_matrix[(MIN (row[i],row[j] )) *gnum+ (MAX (row[i],row[j] ))] -= weight_value;
        }
    }
    
}

//Add offsites back to norm vector and norm matrix
inline void AddOffsites(FLOAT_T *norm_vec, FLOAT_T *norm_matrix, const FLOAT_T  vec_offsite, const FLOAT_T  matrix_offsite, const FLOAT_T  *row_offsites, const FLOAT_T *diag_offsites, uint16 gnum)
{
    
    for (int i = 0; i < gnum; ++i)
    {
        norm_vec[i] += vec_offsite;                     //offsite for norm vector
        
        norm_matrix[i*gnum+i] *= 2;                     //this is doubling diagonal mentioned above
        norm_matrix[i*gnum+i] -= diag_offsites[i];       //we only want to double offsites, but prior values also doubled, should be changed back
        norm_matrix[i*gnum+i] += matrix_offsite;         //offsite for every cell
        
        for (int j = i + 1; j < gnum; ++j)
        {
            
            norm_matrix[i*gnum+j] += row_offsites[i];     //offsite for ith sample
            norm_matrix[i*gnum+j] += row_offsites[j];      //offsite for jth sample
            norm_matrix[i*gnum+j] += matrix_offsite;        //offsite for every cell
            
            norm_matrix[j*gnum + i] = norm_matrix[i*gnum+j];       //square matrix is symmetric
        }
    }
        
}


//Add offsites back to norm vector and norm matrix
inline void AddFlatOffsites(FLOAT_T *norm_vec, FLOAT_T *norm_matrix, const FLOAT_T  vec_offsite, uint16 gnum)
{
    for (int i = 0; i < gnum; ++i)
    {
        norm_vec[i] += vec_offsite;  
        norm_matrix[i*gnum+i] *= 2;                     //this is doubling diagonal mentioned above
        
        for (int j = i + 1; j < gnum; ++j)
        {
            
            norm_matrix[j*gnum + i] = norm_matrix[i*gnum+j];       //square matrix is symmetric
        }
    }
        
}

void flatNorm(const node* tree, const uint16 nodenum, FLOAT_T *norm_vec, FLOAT_T *norm_matrix, const uint16 matrix_size)
{
    vector<uint16> repeat_keys(nodenum, 0 );
    
    uint leaveindex = 0;
    for (int i =0 ; i < nodenum; ++i)
    {
        auto& node = tree[i];
        if (node.numchildren)
        {
            size_t index0 = node.children[0] - tree;
            size_t index1 = node.children[1] - tree;
            
            norm_vec[index0] +=  norm_vec[i];
            norm_vec[index1] +=  norm_vec[i];
            
            for (int j = 0; j < nodenum ; ++j)
            {
                norm_matrix[nodenum*index0+j] +=  norm_matrix[nodenum*i+j];
                norm_matrix[nodenum*index1+j] +=  norm_matrix[nodenum*i+j];
            }
        }
        else
        {
            repeat_keys[leaveindex++] = i;
        }
    }
    
    for (int i =0 ; i < nodenum; ++i)
    {
        auto& node = tree[i];
        if (node.numchildren)
        {
            size_t index0 = node.children[0] - tree;
            size_t index1 = node.children[1] - tree;
            
            
            for (int j = 0; j < nodenum ; ++j)
            {
                norm_matrix[index0+nodenum*j] +=  norm_matrix[i+nodenum*j];
                norm_matrix[index1+nodenum*j] +=  norm_matrix[i+nodenum*j];
            }
        }
    }
    
    for (int i =0 ; i < leaveindex; ++i)
    {
        norm_vec[i] =  norm_vec[repeat_keys[i]];
        
        for (int j = 0 ; j < leaveindex ; ++j)
        {
            norm_matrix[matrix_size*i+j] =  norm_matrix[nodenum*repeat_keys[i]+repeat_keys[j]];
        }
    }
}

/*
void KmerMatrix::getNorm(const uint16* kmervec, const uint16* kmermatrix,const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T  &total_lambda)
{
    //memset(norm_matrix , 0,sizeof(FLOAT_T)*gnum*gnum);
    
    for (size_t i = 0; i < gnum; ++i)
    {
        row_offsites.get()[i] = 0;
        //diag_offsites.get()[i] = norm_matrix[i * gnum + i];
        diag_offsites.get()[i] = 0;
    }
    
    FLOAT_T matrix_offsite = 0, vec_offsite = 0;
    
    FLOAT_T norm_value = 0.0, weight_value = 0.0, mean_repeat = 1.0, mean_repeat_this = 1.0;
    
    uint16 uniqcounter=0, j =0;
    
    float weight = 1.0;
    uint16 lastsize = 0;
    uint16 lastsign = 0;
    const uint16* lastnorm = NULL;
    const uint16* rowdata = kmermatrix;
    
    for (size_t i = 0; i < knum ; ++i)
    {
        
        switch (rowdata[0])
        {
            case '_':
                
                getEachRowValue(depth, kmervec[i], 0, 2,  rowdata[5],1.0 ,total_lambda, norm_value, weight_value);
                break;
            case '=':
                
                mean_repeat = ( (float) rowdata[6] )/255 ;
                mean_repeat *= mean_repeat;
                
                getEachRowValue(depth, kmervec[i], 1, lastsize, rowdata[5], mean_repeat, total_lambda, norm_value, weight_value);
                break;
                
            case '-':
                
                if (lastsign == '-')
                {
                    getEachRowNorm_major(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value, vec_offsite, matrix_offsite, row_offsites.get());
                }
                else if (lastsign == '+')
                {
                    getEachRowNorm_minor(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value, mean_repeat, repeat_keys.get(), repeat_counts.get());
                }
                
                lastsize = rowdata[1]- rowdata[9];
                lastnorm = &rowdata[FIXCOL]+rowdata[9];
                lastsign = rowdata[0];
                norm_value = 0.0;
                weight_value = 0.0;
                
                getEachRowValue(depth, kmervec[i], 0,  2, rowdata[5], 1.0, total_lambda, norm_value, weight_value);
                break;
            case'+':
                
                if (lastsign == '-')
                {
                    getEachRowNorm_major(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value,vec_offsite, matrix_offsite, row_offsites.get());
                }
                else if (lastsign == '+')
                {
                    
                    getEachRowNorm_minor(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value, mean_repeat, repeat_keys.get(), repeat_counts.get());
                    
                }
                
                lastsize = rowdata[1] - rowdata[9];
                lastnorm = &rowdata[FIXCOL+rowdata[9]];
                lastsign = rowdata[0];
                norm_value = 0.0;
                weight_value = 0.0;
                mean_repeat = ( (float) rowdata[6] )/255 ;
                mean_repeat *= mean_repeat;
                
                getEachRowValue(depth, kmervec[i], 1, lastsize, rowdata[5], mean_repeat, total_lambda, norm_value, weight_value);
                break;
            default:
                break;
        }
                    
        rowdata = &rowdata[rowdata[1] + FIXCOL];
        
    }
    
    if (lastsign == '-')
    {
        getEachRowNorm_major(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value,
                             vec_offsite, matrix_offsite, row_offsites.get());
    }
    else if (lastsign == '+')
    {
        getEachRowNorm_minor(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value, mean_repeat, repeat_keys.get(), repeat_counts.get());
    }
    
    AddOffsites(norm_vec, norm_matrix, vec_offsite, matrix_offsite, row_offsites.get(), diag_offsites.get(), gnum);
    
}
*/

void flatPriorNorm2(const node* tree, const uint16 nodenum, FLOAT_T *norm_matrix, const uint16 matrix_size)
{
    vector<uint16> repeat_keys(nodenum, 0 );
    
    uint leaveindex = 0;
    for (int i =0 ; i < nodenum; ++i)
    {
        auto& node = tree[i];
        if (node.numchildren)
        {
            size_t index0 = node.children[0] - tree;
            size_t index1 = node.children[1] - tree;
            
            for (int j = 0; j < nodenum ; ++j)
            {
                norm_matrix[nodenum*index0+j] +=  norm_matrix[nodenum*i+j];
                norm_matrix[nodenum*index1+j] +=  norm_matrix[nodenum*i+j];
            }
        }
        else
        {
            repeat_keys[leaveindex++] = i;
        }
    }
    
    for (int i =0 ; i < nodenum; ++i)
    {
        auto& node = tree[i];
        if (node.numchildren)
        {
            size_t index0 = node.children[0] - tree;
            size_t index1 = node.children[1] - tree;
            
            
            for (int j = 0; j < nodenum ; ++j)
            {
                norm_matrix[index0+nodenum*j] +=  norm_matrix[i+nodenum*j];
                norm_matrix[index1+nodenum*j] +=  norm_matrix[i+nodenum*j];
            }
        }
    }
    
    for (int i =0 ; i < leaveindex; ++i)
    {
        for (int j = 0 ; j < leaveindex ; ++j)
        {
            norm_matrix[matrix_size*i+j] =  norm_matrix[nodenum*repeat_keys[i]+repeat_keys[j]];
        }
    }
}



void KmerMatrix::getNormflat(const uint16* kmervec, const uint16* kmermatrix,const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T  &total_lambda, const node* tree , const uint16 nodenum)
{
    //memset(norm_matrix , 0,sizeof(FLOAT_T)*nodenum*nodenum);
    
    FLOAT_T matrix_offsite = 0, vec_offsite = 0;
    
    FLOAT_T norm_value = 0.0, weight_value = 0.0, mean_repeat = 1.0, mean_repeat_this = 1.0;
    
    uint16 uniqcounter=0, j =0;
    
    uint16 lastsize = 0;
    uint16 lastcount = 0;
    uint16 lastsign = 0;
    uint16 lasttreetotal = 0;
    uint16 lasttreeposi = 0;
    const uint16* lastnorm = NULL;
    const uint16* rowdata = kmermatrix;
    
    float totalchange = 0;
    for (size_t i = 0; i < knum ; ++i)
    {
        switch (rowdata[0])
        {
            case '_':
                getEachRowValue(depth, kmervec[i], 0, 2,  rowdata[5],1.0 ,total_lambda, norm_value, weight_value);
                break;
            case '=':
                getEachRowValue(depth, kmervec[i], 1, lastsize, rowdata[5], mean_repeat, total_lambda, norm_value, weight_value);
                break;
                
            case '-':
               
                getEachRowNorm(lasttreeposi, lasttreetotal,  lastnorm , nodenum, norm_vec,  norm_matrix ,norm_value, weight_value);
                totalchange += weight_value;
                
                lasttreetotal = rowdata[9];
                lasttreeposi = rowdata[8];
                lastsize = rowdata[7];
                lastnorm = &rowdata[FIXCOL];
                lastsign = rowdata[0];
                norm_value = 0.0;
                weight_value = 0.0;
                mean_repeat = 1.0;
                
                getEachRowValue(depth, kmervec[i], 0,  2, rowdata[5], mean_repeat, total_lambda, norm_value, weight_value);
                break;
            case'+':
                
                totalchange += weight_value;
                getEachRowNorm(lasttreeposi, lasttreetotal,  lastnorm , nodenum, norm_vec,  norm_matrix ,norm_value, weight_value);
                
                lasttreetotal = rowdata[9];
                lasttreeposi = rowdata[8];
                lastsize = rowdata[7];
                lastnorm = &rowdata[FIXCOL];
                lastsign = rowdata[0];
                norm_value = 0.0;
                weight_value = 0.0;
                mean_repeat = ( (float) rowdata[6] )/255 ;
                mean_repeat *= mean_repeat;
                
                getEachRowValue(depth, kmervec[i], 1, lastsize, rowdata[5], mean_repeat, total_lambda, norm_value, weight_value);
                break;
            default:
                break;
        }
                    
        rowdata = &rowdata[rowdata[1] + FIXCOL];
        
    }
    
    totalchange += weight_value;
    
    
    getEachRowNorm(lasttreeposi, lasttreetotal, lastnorm, nodenum, norm_vec,  norm_matrix,norm_value, weight_value);
    AddFlatOffsites(norm_vec, norm_matrix, vec_offsite, nodenum);
    flatNorm(tree,nodenum, norm_vec, norm_matrix, gnum);
}






