//
//  Genotyper.hpp
//  Ctyper2
//
//  Created by walfred on 5/20/25.
//

#ifndef Genotyper_hpp
#define Genotyper_hpp

#include <stdio.h>
#include <vector>
#include <utility>
#include <unordered_set>

#include "config.hpp"
#include "KmerCounter.hpp"
#include "Regression.hpp"
#include "TreeRound.hpp"
#include "KtableReader.hpp"
#include "PriorData.hpp"
#include "KmerWindow.hpp"
#include "KmerMatrix.hpp"

#define DefaultSize 5000
#define DefaultKmeralloc 1000000


class Genotyper
{
    using kmer_int = ull;
    
public:
    
    unique_ptr<int[]> results;
    unique_ptr<FLOAT_T[]> reminders;
    unique_ptr<FLOAT_T[]> coefs;
    unique_ptr<FLOAT_T[]> residuels;
    const size_t knum, pnum;
    const uint window;
    const int Nsubthreads;
    Genotyper(size_t k, size_t p, KmerCounter &c ,PriorData &priordata, const int w, const int N):
    knum(k),
    pnum(p),
    window(w),
    counter(c),
    priordata_manager(priordata),
    all_kmer_counts(new counterint[k+1]),
    kmer_counts(new uint16[DefaultKmeralloc]),
    Nsubthreads(N),
    
    norm_vec(new FLOAT_T[MEMEXP*DefaultSize]),
    norm_matrix(new FLOAT_T[MEMEXP*MEMEXP *DefaultSize*DefaultSize]),
    reduce_matrix(new FLOAT_T[DefaultSize*DefaultSize]),
    
    coefs(new FLOAT_T[MAX_UINT16]),
    residuels(new FLOAT_T[MAX_UINT16]),
    reminders(new FLOAT_T[MAX_UINT16]),
    results(new int[MAX_UINT16]),
    
    finished_group(p)
    {};
    
    void counting(const std::string& inputfile)
    {

        cerr << "counting kmers for sample: " << inputfile<<endl;

        auto begin = std::chrono::high_resolution_clock::now();
                
        counter.Call(inputfile.c_str(), all_kmer_counts.get(), excess_kmers, totalbases, totalreads, totalbgs, Nsubthreads);

        finishcounting = 1;

        auto end = std::chrono::high_resolution_clock::now();
        
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            
        cerr<<"finished counting "<< inputfile <<" at time: "<<elapsed.count()* 1e-9 <<endl;
        
    };
    
    void transfercounts(const vector<uint>& hashs, const size_t numhash)
    {
        for (int i = 0; i < numhash; ++i)
        {
            auto index = hashs[i];
            uint16 count = all_kmer_counts[index];
            
            if ( __builtin_expect(count < 255, 1))
            {
                kmer_counts.get()[i] = count;
            }
            else
            {
                kmer_counts.get()[i] = count + excess_kmers.findvalue(hashs[i]);
            }
        }
    }
    
    void runOneGroup(const PriorChunk* priorData, const std::string& inputfile, const std::string& outputfile, const float depth, std::mutex& Threads_lock)
    {
        auto t1 = std::chrono::steady_clock::now();
        
        transfercounts(priorData->kmerhashs, priorData->kmervec_size);
        
        cout << "Starting " << priorData->prefix << " for sample: " << inputfile <<endl;
        matrix.getNormflat(kmer_counts.get(), priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size, norm_vec.get(), norm_matrix.get(), total_lambda,  priorData->phylo_tree, priorData->nodenum);
        //matrix.getNorm(kmer_counts.get(), priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size, norm_vec.get(), norm_matrix.get(), total_lambda);
        
        
        
        auto t2 = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Finished generating matrix for sample: " << inputfile
                  << " at " << elapsed << " ms" << std::endl;
        
        regresser.Call(kmer_counts.get(), priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size, norm_vec.get(), norm_matrix.get(), reduce_matrix.get(), total_lambda, priorData->gene_kmercounts, coefs.get(), residuels.get(), priorData->numgroups, priorData->genegroups, priorData->numsmallgroups, priorData->smallgroups, priorData->groupkmernums);

        
        auto t3 = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
        std::cout << "Finished regerssion for sample: " << inputfile
                  << " at " << elapsed << " ms" << std::endl;
        
        cout << "rounding for sample: " << inputfile<<endl;
        
        tree.Run(priorData->phylo_tree, coefs.get(), gnum, &results.get()[0], &reminders.get()[0], residuels.get(), norm_matrix.get());

        cout << "determine window residuels: " << inputfile<<endl;

    
        KmerWindow kmerwindow(window);
        //kmerwindow.resize(priorData->pathsizes);
        //kmerwindow.WindowCovers(kmer_counts.get(), priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size, priorData->genenum, &results.get()[0], total_obs, total_exp);
        
        vector<vector<tuple<int, int, float, string>>> PatialCopies(priorData->pathnames.size()+1);
        //kmerwindow.PartialCopy(PatialCopies, &reminders.get()[0], priorData->genenames, priorData->pathnames, depth);
        
        write(priorData, outputfile, inputfile, priorData->prefix, priorData->genenames, PatialCopies, kmerwindow.windowcovers, depth, Threads_lock);

        cout<<"finish run"<<endl;
    };
    
    void write(const PriorChunk* priorData, const std::string& outputfile, const string &sample, const string &prefix, const vector<string>&genenames_ori, const vector<vector<tuple<int, int, float, string>>>& PatialCopies, const vector<vector<tuple<int,int,int>>>& windowcovers, const float depth, std::mutex& Threads_lock)
    {
        
        auto genenames(genenames_ori);
        for (auto &genename: genenames)
        {
            genename = genename.substr(0,genename.find('\t', 0));
        }
        
        std::unique_lock<std::mutex> lck(Threads_lock);
        
        FILE *fwrite;
        
        if (outputfile != "stdout")
        {
            fwrite=fopen(outputfile.c_str(), "a");
        }
        else
        {
            fwrite=stdout;
        }
        
        if (fwrite==NULL)
        {
            std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
            
            std::_Exit(EXIT_FAILURE);
        }
        
        fprintf(fwrite,">%s\t%s\n", prefix.c_str(), sample.c_str());
        //fprintf(fwrite,"rsdl: %.4lf\n", regress.first);
        fprintf(fwrite,"lambda: %llu/%llu\n",total_obs, total_exp);
        
        fprintf(fwrite,"regress: ");
        const float cutoff = 0.5 / (gnum + 1);
        for (int i = 0; i < gnum; ++i)
        {
            if (coefs.get()[i] > cutoff) fprintf(fwrite,"%s:%.2lf,", genenames[i].c_str(),coefs.get()[i]);
        }
        fprintf(fwrite,"\n");
    
        fprintf(fwrite,"round: ");
        for (int i = 0; i < gnum; ++i)
        {
            int result = results.get()[i];
            
            if (result > 0)
            {
                
                for (int j = 0 ; j < result  ; ++j)
                {
                    string info = "";
                    if (priorData->gene_kmercounts[i] < 1000) info = "(aux)";
        
                    fprintf(fwrite,"%s%s,", genenames[i].c_str(), info.c_str());
                }
            }
        }
        fprintf(fwrite,"\n");
        
        fprintf(fwrite,"result: ");
        for (int i = 0; i < gnum; ++i)
        {
            int result = results.get()[i];
            
            if (result  > 0)
            {
                for (int j = 0 ; j < result ; ++j)
                {
                    if (priorData->gene_kmercounts[i] >= 1000) fprintf(fwrite,"%s,", genenames[i].c_str());
                }
            }
        }
        fprintf(fwrite,"\n");
        
        fclose(fwrite);
        
        return ;
    };

    void newsample()
    {
        memset(all_kmer_counts.get(), 0, sizeof(counterint) * knum);
        excess_kmers.reset();
        finished_group.assign(pnum , 0);
        finishcounting = 0;
        totalbases = 0;
        totalreads = 0;
        totalbgs = 0;
    }
    
    void newgroup(const PriorChunk* priorData)
    {
        if (MAX(gnum + 10 , alloc_size) > DefaultSize)
        {
            alloc_size = MAX( DefaultSize, gnum + 10 );
            
            norm_vec.reset(new FLOAT_T[MEMEXP*alloc_size]);
            
            try_allocate_unique(norm_matrix, MEMEXP*MEMEXP*alloc_size*alloc_size, MEMEXP*MEMEXP*alloc_size*alloc_size);
                        
            coefs.reset(new FLOAT_T[alloc_size]);
            
            residuels.reset(new FLOAT_T[alloc_size]);
            
            results.reset(new int[alloc_size]);
        }

        if (MAX( priorData->numsmallgroups + 10 , alloc_size2) > DefaultSize)
        {
            alloc_size2 = MAX( DefaultSize, priorData->numsmallgroups + 10 );
            
            try_allocate_unique(reduce_matrix, alloc_size2*alloc_size2, alloc_size2*alloc_size2);

        }
        
        if (MAX(priorData->kmervec_size + 10 , kmer_alloc_size) > DefaultKmeralloc)
        {
            kmer_alloc_size = MAX ( priorData->kmervec_size + 10 , DefaultKmeralloc);

            kmer_counts.reset(new uint16[kmer_alloc_size]);
        }
        
        memset(kmer_counts.get(), 0, sizeof(uint16) * priorData->kmervec_size);
        
        memcpy(norm_matrix.get(), priorData->prior_norm, sizeof (FLOAT_T) *  priorData->prior_norm_allocsize);
                                
        memset(norm_vec.get(), 0, sizeof (FLOAT_T) * MEMEXP * gnum );
        
        memset(coefs.get(), 0, sizeof(FLOAT_T) * gnum);
        
        memset(residuels.get(), 0, sizeof(FLOAT_T) * gnum);
        
        memset(results.get(), 0, sizeof(int) * gnum);
        
        //kmerwindow.resize(priorData->pathsizes);
                        
        total_lambda = 0;
        total_obs = 0;
        total_exp = 0;
                
    };
    
    void run(const std::string& inputfile, const std::string& outputfile, float depth,std::mutex& Threads_lock)
    {
        cerr<<"running for sample: "<<inputfile << endl;

        newsample();
        
        counting(inputfile);
        
        if (depth <= 0)
        {
            depth = ( 0.5 * totalbgs )/counter.totalbks;
        }
        
        FILE *fwrite;
        
        if (outputfile != "stdout")
        {
            fwrite=fopen(outputfile.c_str(), "a");
        }
        else
        {
            fwrite=stdout;
        }
        
        if (fwrite==NULL)
        {
            std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
            std::_Exit(EXIT_FAILURE);
        }
                
        fprintf(fwrite,"@totalreads: %llu, totalbackgrounds: %llu/%llu \n", totalreads, totalbgs, counter.totalbks);
        
        fclose(fwrite);
        
        for (int i = 0; i < pnum; ++i)
        {
            
            auto begin = std::chrono::high_resolution_clock::now();
            
            PriorChunk* priorData = priordata_manager.getNextChunk(finished_group);
            
            cout << "running gene " << priorData->prefix << " for sample " << inputfile << endl ;
            
            gnum = priorData->genenum;
            
            newgroup(priorData);
                        
            runOneGroup (priorData, inputfile, outputfile, depth, Threads_lock);
    
            finished_group[priorData->index] = 1;
            
            priordata_manager.FinishChunk(priorData);

            auto end = std::chrono::high_resolution_clock::now();

            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            
            cout<<"finished sample:" << inputfile << "for gene:"<< priorData->prefix << " for :" << elapsed.count()* 1e-9 <<endl;
        }

    cerr<<"finished sample: "<<inputfile << endl;
    };
    
    

private:
    
    unique_ptr<uint16[]> kmer_counts;
    unique_ptr<counterint[]> all_kmer_counts;
    Excess_hash excess_kmers;
    unique_ptr<FLOAT_T[]> norm_vec;
    unique_ptr<FLOAT_T[]> norm_matrix;
    unique_ptr<FLOAT_T[]> reduce_matrix;
    
    FLOAT_T total_lambda =0;
    ull total_exp = 0, total_obs = 0;
    
    PriorData &priordata_manager;
    KmerCounter &counter;
    KmerMatrix matrix;
    Regression regresser;
    TreeRound tree;
    //KmerWindow kmerwindow;
    
    ull totalbases = 0, totalreads = 0, totalbgs = 0;
    size_t gnum;
    size_t alloc_size = DefaultSize;
    size_t alloc_size2 = DefaultSize;
    size_t kmer_alloc_size = DefaultKmeralloc;
    bool finishcounting = 0;
    vector<bool> finished_group;
    
};


#endif /* Genotyper_hpp */



