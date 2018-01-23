//
//  main.cpp
//  quotientFilter
//
//  Created by Matthew Green on 8/30/17.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include "hashType.h"

vector<string> getKmersFromString(string read, int window_size)
{
    vector<string> toRet;
    string kmer = "";
    
    for(int j = 0; j < window_size; j++)
    {
        kmer += read[j];
    }
    
    toRet.push_back(kmer);
    
    for(int i = 1; i < read.size() - window_size + 1; i++)
    {
        
        kmer = kmer.substr(1);
        kmer += read[i + window_size - 1];
        
        toRet.push_back(kmer);
    }
    return toRet;
}



int main(int argc, char *argv[])
{
    if(argc < 4)
    {
        cout << "Place arguments in the order of 1. Q size, 2. R size, 3. Fasta File" << endl;
    }
    //recommended q = 28, R = 16
    quotientFilter *filter = new quotientFilter(atoi(argv[1]), atoi(argv[2]));
    ifstream fastaFile(argv[3]);
    string line;
    char kmer[WIN_SIZE_PLUS];
    vector<string> lineKmers;
    int i = 0;
    double numKmers = 0;
    
    unordered_map<string, int> kmers;
    
    while(getline(fastaFile, line))
    {
        if(i % 2 == 1)
        {
            lineKmers.clear();
            lineKmers = getKmersFromString(line,32);
            numKmers += lineKmers.size();
            for(auto k : lineKmers)
            {
                //kmers[k] = 1;
                strcpy(kmer, k.c_str());
                filter->insertKmer(kmer);
            }
        }
        if(i % 5000 == 0)
        {
            cout << "Processed: " << i/2 << endl;
        }
        
        i++;
        
    }
    
    cout << "Loaded: " << numKmers << " Kmers" << endl;
    
    filter->build_S();
    
    
    fastaFile.close();
    
    return 0;
}
