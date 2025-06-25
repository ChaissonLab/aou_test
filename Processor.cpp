
#include "Processor.hpp"

void Processor::Run()
{
    
    if (genes.size() > 0)
    {
        totalgroups = priordata_manager.LoadIndex(genes);
    }
    else
    {
        totalgroups = priordata_manager.LoadIndex();
    }
   
    cout<<"reading all kmer targets"<<endl;
    
    //kmer_hash.initiate(4 * priordata_manager.totalkmers);
 
    Counter = new KmerCounter;
    if (backgroundfile.length() > 0) Counter->load_backgrounds(backgroundfile.c_str());
        
    //Counter->LoadRegion(regions);
    totalkmers = Counter->read_target(matrixfile.c_str());
    
    
    priordata_manager.Addhash(&(Counter->kmer_hash));
    cout<<"finishing reading targets, start genotyping"<<endl;

    std::vector<std::unique_ptr<std::thread>> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        threads.push_back(std::unique_ptr<std::thread>(new std::thread(&Processor::Onethread, this)));
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i].get()->join();
    }
    
    
    
}

void Processor::Onethread()
{
    
    unique_ptr<Genotyper> genotyper = unique_ptr<Genotyper>(new Genotyper(totalkmers, totalgroups,  *Counter, priordata_manager, window, Nsubthreads));
    
    
    while (restfileindex < inputfiles.size() )
    {
        while(Threads_lock.try_lock());
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size() ) break;
        
    
        string outputfile;
        if (inputindex < outputfiles.size())
        {
            outputfile = outputfiles[inputindex];
        }
        else if (outputfiles.size())
        {
            outputfile = outputfiles[outputfiles.size()-1];
        }
        else
        {
            outputfile = "stdout";
        }
                
        float depth = -1;
        if (inputindex < depths.size())
        {
            depth= depths[inputindex];
        }
        else if (depths.size() > 0)
        {
            depth= depths[depths.size() -1];
        }
        
        genotyper.get()->run(inputfiles[inputindex], outputfile, depth, Threads_lock2);
    }
    
    
}

