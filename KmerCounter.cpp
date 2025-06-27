
#include "KmerCounter.hpp"

extern bool ifbg;

template <typename T>
static void kmer_deconpress(const char * kmer_compress, T &larger_kmer)
{
    larger_kmer = 0;
    for (int pos = 0; pos < 30; ++pos)
    {
        if (kmer_compress[pos] == '\t') break;
        larger_kmer <<= 6;
        larger_kmer += kmer_compress[pos] - '0';
    }
    
}

template <typename T>
inline static void update_counter(Kmer_hash &Kmer_hash,ull &larger_kmer, T* vec, Excess_hash &exbucks)
{
    
    Kmer_hash.addvalue(larger_kmer, vec, exbucks);

	/*
    if (index != NULL)
    {
        if ( *index == 0) return;
        
        if ( __builtin_expect(vec[*index] < 255, 1))
        {
            vec[*index] ++;
        }
        else
        {
            exbucks.add(*index, 1);
        }
    }
	*/
}




template <typename T>
static void kmer_read_31(char base,  std::size_t &current_size, T &current_kmer, T &reverse_kmer)
{
    int converted = 0;
    T reverse_converted;
    const int klen = 31;
    
    const T mask = ((T)1 << (klen*2)) - 1;
    const int shift = 2*klen - 2;
    
    //if (base == '\n' || base == ' ') return;
    
    if (base_to_int(base, converted))
    {
        
        current_kmer = ((current_kmer << 2) & mask) | converted;
        reverse_kmer = (reverse_kmer >> 2) | ((T)(0b11 ^ converted) << shift);
        
    }
    
    else
    {
        current_size = -1;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

/*
template <int dictsize>
void KmerCounter<dictsize>::LoadRegion(std::vector<char *> &r)
{
    this->regions = r;
}
*/

void KmerCounter::load_backgrounds(const char * backgroudfile)
{
    FastaReader fastafile(backgroudfile);
    
    backgrounds.resize(bgkmernum,0);
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
        
    std::hash<std::uint64_t> hasher;
    std::string StrLine;
    StrLine.resize(MAX_LINE);
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
            
            if (base == '\n' || base == ' ') continue;

            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen || base >= 'a') continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
                        
            uint hash_ = larger_kmer % bgkmernum;
			
			if (hash_ != 0)
			{
				backgrounds[hash_] = larger_kmer;
				totalbks ++;
			}
        }
    }
    
    fastafile.Close();
    
}

ull KmerCounter::read_target(FastaReader &fastafile)
{
    
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
        
    std::string StrLine;
    StrLine.resize(MAX_LINE);
   
    uint hash_;
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
            
            if (base == '\n' || base == ' ') continue;

            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen || base >= 'a') continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            kmer_hash.preview(larger_kmer);
                
        }
    }
    
    fastafile.Reset();
    kmer_hash.initiate();
    
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
            
            if (base == '\n' || base == ' ') continue;

            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen || base >= 'a') continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            totalkmers += kmer_hash.add(larger_kmer);
                
        }
    }
    
    return kmer_hash.totaleles + 10;

};

ull KmerCounter::read_target(KtableReader &ktablefile)
{
    
    kmer_int current_kmer = 0;
        
    std::string StrLine;
    StrLine.resize(MAX_LINE);
    
    char base;
    int pos;
    int elecounter;
    while (ktablefile.nextLine_kmer(StrLine))
    {
        if (StrLine[0] == '^') continue;
        pos = 0;
        
        elecounter = 0;
        for (; pos < klen; ++pos)
        {
            base = StrLine[pos];
            if (base == '\t' && ++elecounter > FIXCOL-3) break;
        }
        ++pos;
                
        kmer_deconpress(&StrLine.c_str()[pos], current_kmer);
        
        kmer_hash.preview(current_kmer);
        
    }
    
    ktablefile.Reset();
    
    kmer_hash.initiate();
	int i = 0;
    while (ktablefile.nextLine_kmer(StrLine))
    {
        if (StrLine[0] == '^') continue;
        pos = 0;
        elecounter = 0;
        for (; pos < klen; ++pos)
        {
            base = StrLine[pos];
            if (base == '\t' && ++elecounter > FIXCOL-3) break;
        }
        ++pos;
                
        kmer_deconpress(&StrLine.c_str()[pos], current_kmer);
        
        totalkmers += kmer_hash.add(current_kmer);
        
    }
    
    return kmer_hash.totaleles + 10;
        
};
 
ull KmerCounter::read_target(const char* inputfile)
{
    int pathlen = (int)strlen(inputfile);
    
    size_t totalhash = 0;
    if ( (pathlen > 2 && strcmp(inputfile+(pathlen-3),".fa")==0) || (pathlen > 6 && strcmp(inputfile+(pathlen-6),".fasta") == 0 ))
    {
        FastaReader readsfile(inputfile);
        totalhash = read_target(readsfile);
    }
    else
    {
        
        KtableReader readsfile(inputfile);
        totalhash = read_target(readsfile);
    }
    
    return totalhash;
}

void KmerCounter::Call(const char* inputfile, counterint* samplevecs, Excess_hash &exbucks, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads, std::string reference)
{
    
    vector<string> files;
    if (std::filesystem::is_directory(inputfile))
    {
        for (const auto& entry : std::filesystem::directory_iterator(inputfile))
        {
            if (std::filesystem::is_regular_file(entry.path()))
            {
                files.push_back(entry.path().string());
            }
        }
    }
    else
    {
        files.push_back(inputfile);
    }
    
    for (string& filestring: files)
    {
        const char* file = filestring.c_str();
        
        int pathlen = (int)strlen(file);
        
        if ( pathlen > 2 && strcmp(file+(pathlen-3),".gz") == 0 )
        {
            
            FastqReader readsfile(file);
            count_kmer_(readsfile, samplevecs, exbucks, nBases, nReads, nBg, nthreads);
            readsfile.Close();
        }
            
        else if ( (pathlen > 2 && strcmp(file+(pathlen-3),".fa")==0) || (pathlen > 6 && strcmp(file+(pathlen-6),".fasta") == 0 ))
        {
            FastaReader readsfile(file);
            
            count_kmer_(readsfile, samplevecs, exbucks, nBases, nReads, nBg, 1);
            readsfile.Close();
        }
        else if (pathlen > 5 && ( strcmp(file+(pathlen-5),".cram") == 0  || strcmp(file+(pathlen-4),".bam") == 0 || strcmp(file+(pathlen-4),".sam") == 0  ))
        {
	  CramReader readsfile(file, reference);
            readsfile.LoadRegion(regions);
            count_kmer_(readsfile, samplevecs, exbucks, nBases, nReads, nBg, nthreads);
            readsfile.Close();
        }
    }
    
}

template <class typefile>
void KmerCounter::count_kmer_(typefile &file, counterint* samplevecs, Excess_hash &exbucks, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads)
{
    
    std::vector<std::thread> threads;
        
    for(int i=0; i< nthreads; ++i)
    {
        threads.push_back(std::thread([this, &file, &samplevecs, &exbucks, &nBases, &nReads, &nBg]()
        {
            this->count_kmer(file, samplevecs,exbucks, nBases, nReads, nBg);
        }));
    }
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i].join();
    }
    
};



void KmerCounter::count_kmer(FastaReader &file, counterint* samplevecs, Excess_hash &exbucks, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    uint hash_;
    std::string StrLine;
    std::string Header;
    
    while (file.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                current_kmer = 0;
                reverse_kmer = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        for (auto base: StrLine)
        {
            
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ') continue;
 
            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            update_counter(kmer_hash, larger_kmer, samplevecs, exbucks);
                        
            if ( ifbg && __builtin_expect(backgrounds[larger_kmer % bgkmernum] == larger_kmer , 0)  )
            {
                nBg ++;
            }
            
        }
        
    //nBases+=MAX(0, StrLine.length() - klen + 1);
    nReads+=1;
    if (nReads % 10000000 == 0) {
      cerr << "processed " << nReads / 1000000 << "M reads." << endl;
    }
    }
        
    
    return ;
};


void KmerCounter::count_kmer(CramReader &file, counterint* samplevecs, Excess_hash &exbucks, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string Header;
    
    uint hash_ ;
    size_t rlen = 0;
    uint8_t *StrLine = NULL;
    bool ifbackground = 0;
    char base;
    while (file.nextLine_prt(StrLine, rlen))
    {
        for (int pos = 0; pos < rlen; ++pos)
        {
        base = seq_nt16_str[bam_seqi(StrLine,pos)];
            //if (base=='\n') break;
            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;

            update_counter(kmer_hash, larger_kmer, samplevecs, exbucks);
           
            if ( ifbg && __builtin_expect(backgrounds[larger_kmer % bgkmernum] == larger_kmer , 0)  )
            {
                nBg ++;
            }

        }
        
        //nBases+=MAX(0, StrLine.length() - klen + 1);
        nReads+=1;
        if (nReads % 10000000 == 0){
            cerr << "read " << nReads / 1000000 << "M reads." << endl;
        }
    }
    cerr << "Finished reading " << nReads << " reads." << endl; 
    return ;
};

void KmerCounter::count_kmer(FastqReader &file, counterint* samplevecs, Excess_hash &exbucks, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string Header;
    
    uint hash_ = 0;
    size_t rlen = 0;
    const char *StrLine = NULL;
    
    while (file.nextLine_prt(StrLine, rlen))
    {
        for (int pos = 0; pos < rlen; ++pos)
        {
            char base = StrLine[pos];
            
            //if (base=='\n') break;
            
            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            update_counter(kmer_hash,  larger_kmer, samplevecs, exbucks);
            
            if ( ifbg && __builtin_expect(backgrounds[larger_kmer % bgkmernum] == larger_kmer , 0)  )
            {
                nBg ++;
            }
            
            //nBases+=MAX(0, StrLine.length() - klen + 1);
            nReads+=1;
            if (nReads % 10000000 == 0) {
                cerr << "processed " << nReads / 1000000 << "M reads." << endl;
            }
        }
        
        return ;
    }
}


