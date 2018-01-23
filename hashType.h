//
//  hashType.h
//  quotientFilter
//
//  Created by Matthew Green on 8/30/17.
//

#ifndef quotientFilter_hashType_h
#define quotientFilter_hashType_h

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <string>
#include <queue>
#include <cstdlib>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <ctime>

#define BP_PER_UCHAR 4
#define BITS_PER_UCHAR 8
#define WIN_SIZE 32
#define K_SIZE (WIN_SIZE+1)/BP_PER_UCHAR
#define WIN_SIZE_PLUS 33
#define MOD 2147483647
#define HL 31
#define IS_OCCUPIED 0x4
#define IS_CONTINUATION 0x2
#define IS_SHIFTED 0x1
#define SET_OCCUPIED_FALSE 0x3
#define SET_CONTINUATION_FALSE 0x5
#define SET_SHIFTED_FALSE 0x6
#define USED 0x80

using namespace std;

struct decodeElement
{
    unsigned char remainder[2];
    unsigned short count;
    unsigned char used;
}decodeElement;

struct QF_Entry
{
    unsigned char bitHolder;
    
    unsigned char remainder[2];
    unsigned short count;
}QF_Entry;

struct decodeInfo
{
    int decodeSize;
    int clusterStartLocation;
    int clusterSize;
}decodeInfo;

class quotientFilter
{
private:
    vector<vector<struct decodeElement>> decoded;
    int q;
    int r;
    long size;
    struct QF_Entry *tab_entry;
    
    void initDecodedStructure()
    {
        vector<struct decodeElement> empty;
        for(int i = 0; i < 25000; i++)
        {
            decoded.push_back(empty);
        }
        
        cout << "Decoded Structure initiated!" << endl;
    }
    
    void clearDecodedStructure(int toClear)
    {
        for(int i = 0; i < toClear; i++)
        {
            decoded[i].clear();
        }
    }
    
    unsigned createMask(unsigned a, unsigned b)
    {
        unsigned r = 0;
        for (unsigned i=a; i<=b; i++)
            r |= (1 << i);
        
        return r;
    }
    
    long convert_bin_to_hash(unsigned char *result, int q)
    {
        long toRet = 0x0000, extra = 0x0000;
        int temp = 0, iterations = q / 8, extraBits = q % 8;
        for(int i = 0; i < iterations; i++)
        {
            q -= 8;
            temp = 0;
            temp = result[i] << q;
            toRet |= temp;
        }
        
        if(extraBits)
        {
            extra = result[iterations] >> (8 - extraBits);
            
            toRet |= extra;
        }
        
        return toRet;
    }
    
    unsigned char* convert_bin_to_remainder(unsigned char *result, int q, int r)
    {
        int loc = r/8;
        unsigned char *toRet = result + 8 - loc;
        
        return toRet;
    }
    
    void encode(struct decodeInfo info)
    {
        long curLoc = info.clusterStartLocation;
        
        for(int i = 0; i < info.decodeSize; i++)
        {
            for(int j = 0; j < decoded[i].size(); j++)
            {
                //always set the is_occupied flag to true
                this->tab_entry[i + info.clusterStartLocation].bitHolder |= IS_OCCUPIED;
                memcpy(this->tab_entry[curLoc].remainder, decoded[i][j].remainder, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
                this->tab_entry[curLoc].count = decoded[i][j].count;
                
                //we are placing the anchor of the cluter in so there is no shift or continuation
                if(i == 0 && j == 0)
                {
                    this->tab_entry[curLoc].bitHolder &= SET_SHIFTED_FALSE;
                    this->tab_entry[curLoc].bitHolder &= SET_CONTINUATION_FALSE;
                }
                else
                {
                    this->tab_entry[curLoc].bitHolder |= IS_SHIFTED;
                    
                    //if j == 0 then that means that this is the beginning of a run, if j > 0, then
                    //encoded[curLoc] is part of a run, so set the is_continuation flag
                    if(j > 0)
                    {
                        this->tab_entry[curLoc].bitHolder |= IS_CONTINUATION;
                    }
                    else
                    {
                        this->tab_entry[curLoc].bitHolder &= SET_CONTINUATION_FALSE;
                    }
                }
                
                curLoc++;
                
            }
        }
    }
    
    struct decodeInfo decodeCluster(int clusterLoc)
    {
        std::queue<int> canonicalQueue;
        int listLoc = 0, startClusterLoc;
        struct decodeInfo toRet;
        struct decodeElement temp;
        
        //here we find the beginning of the cluster
        while(this->tab_entry[clusterLoc].bitHolder & IS_SHIFTED)
        {
            clusterLoc--;
        }
        
        startClusterLoc = clusterLoc;
        
        //place the "anchor" of the cluster into the data structure
        memcpy(temp.remainder, this->tab_entry[clusterLoc].remainder, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
        temp.count = this->tab_entry[clusterLoc].count;
        temp.used = this->tab_entry[clusterLoc].bitHolder & USED;
        decoded[listLoc].push_back(temp);
        
        //listLoc++;
        clusterLoc++;
        
        //now start adding to the data structure with all of the shifted elements aka members of the cluster
        while(this->tab_entry[clusterLoc].bitHolder & IS_SHIFTED)
        {
            
            //whenever we find an occupied position, add it to the queue
            if(this->tab_entry[clusterLoc].bitHolder & IS_OCCUPIED)
            {
                canonicalQueue.push(clusterLoc);
            }
            
            //if we are not at a continuation then we are at the beginning of a new run, so add
            //to the location being held in our queue
            if(!(this->tab_entry[clusterLoc].bitHolder & IS_CONTINUATION))
            {
                listLoc = canonicalQueue.front();
                canonicalQueue.pop();
                listLoc -= startClusterLoc;
            }
            memcpy(temp.remainder, this->tab_entry[clusterLoc].remainder, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
            temp.count = this->tab_entry[clusterLoc].count;
            temp.used = this->tab_entry[clusterLoc].bitHolder & USED;
            
            //push back our "remainder" to the location being worked on
            decoded[listLoc].push_back(temp);
            
            
            clusterLoc++;
        }
        
        toRet.clusterSize = clusterLoc - startClusterLoc;
        toRet.clusterStartLocation = startClusterLoc;
        toRet.decodeSize = listLoc + 1;
        
        //return the starting location in the cluster
        return toRet;
    }
    
    struct decodeInfo insertDecodeCluster(int clusterLoc)
    {
        std::queue<int> canonicalQueue;
        int listLoc = 0, startClusterLoc;
        struct decodeInfo toRet;
        struct decodeElement temp;
        vector<struct decodeElement> empty;
        
        //here we find the beginning of the cluster
        while(this->tab_entry[clusterLoc].bitHolder & IS_SHIFTED)
        {
            clusterLoc--;
        }
        
        startClusterLoc = clusterLoc;
        
        //place the "anchor" of the cluster into the data structure
        memcpy(temp.remainder, this->tab_entry[clusterLoc].remainder, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
        temp.count = this->tab_entry[clusterLoc].count;
        decoded[listLoc].push_back(temp);
        
        //listLoc++;
        clusterLoc++;
        
        //now start adding to the data structure with all of the shifted elements aka members of the cluster
        while(((this->tab_entry[clusterLoc].bitHolder & IS_SHIFTED) || (this->tab_entry[clusterLoc].bitHolder & IS_OCCUPIED)) && (clusterLoc < this->size))
        {
            
            //whenever we find an occupied position, add it to the queue
            if(this->tab_entry[clusterLoc].bitHolder & IS_OCCUPIED)
            {
                canonicalQueue.push(clusterLoc);
            }
            
            //if we are not at a continuation then we are at the beginning of a new run, so add
            //to the location being held in our queue
            if(!(this->tab_entry[clusterLoc].bitHolder & IS_CONTINUATION))
            {
                listLoc = canonicalQueue.front();
                canonicalQueue.pop();
                listLoc -= startClusterLoc;
            }
            
            memcpy(temp.remainder, this->tab_entry[clusterLoc].remainder, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
            temp.count = this->tab_entry[clusterLoc].count;
            
            while(decoded.size() < listLoc + 2)
            {
                decoded.push_back(empty);
            }
            
            //push back our "remainder" to the location being worked on
            decoded[listLoc].push_back(temp);
            
            clusterLoc++;
        }
        
        toRet.clusterSize = clusterLoc - startClusterLoc;
        toRet.clusterStartLocation = startClusterLoc;
        toRet.decodeSize = listLoc + 2;
        
        //return the starting location in the cluster
        return toRet;
    }
    
    bool compareRemainders(unsigned char *r1, unsigned char *r2, int sizeR)
    {
        for(int i = 0; i < sizeR / 8; i++)
        {
            if(r1[i] ^ r2[i])
            {
                return false;
            }
        }
        
        return true;
    }
    
    void insertElement(unsigned int p, unsigned char *r, unsigned char *result)
    {
        int decodedIndex;
        struct decodeElement temp;
        bool toEncode = true;
        
        //if we are inserting into an empty location, then we can just insert without decoding
        if(!(this->tab_entry[p].bitHolder & IS_OCCUPIED) && !(this->tab_entry[p].bitHolder & IS_SHIFTED))
        {
            this->tab_entry[p].bitHolder |= IS_OCCUPIED;
            memcpy(this->tab_entry[p].remainder, r, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
            this->tab_entry[p].count = 1;
        }
        else
        {
            struct decodeInfo info;
            vector<struct decodeElement> empty;
            //CHANGED: made this into insertDecodeCluster
            info = insertDecodeCluster(p);
            decodedIndex = p - info.clusterStartLocation;
            
            while(decoded.size() < decodedIndex + 1)
            {
                decoded.push_back(empty);
                info.decodeSize = (int)decoded.size();
            }
            
            //if where i need to encode will exceed the QF range, then ignore
            if(info.clusterStartLocation + info.clusterSize >= this->size - 1)
            {
                toEncode = false;
            }
            else if(decoded[decodedIndex].size() == 0)// || decoded[decodedIndex].size() == 1)
            {
                //cout << "zero: " << decodedIndex << " " << info->decodeSize << endl;
                memcpy(temp.remainder, r, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
                //temp.remainder = r;
                temp.count = 1;
                decoded[decodedIndex].push_back(temp);
                //this if statement seems to be the problem of some bugs... perhaps im calculating decode size wrong
                if(decodedIndex + 1 > info.decodeSize)
                {
                    info.decodeSize = decodedIndex + 1;
                }
            }
            else
            {
                vector<struct decodeElement>::iterator it, end;
                int k = 0;
                it = decoded[decodedIndex].begin();
                end = decoded[decodedIndex].end();
                //cout << decoded[decodedIndex].size() << endl;
                
                while(r[0] > decoded[decodedIndex][k].remainder[0] && !compareRemainders(r, decoded[decodedIndex][k].remainder, this->r) && k + 1 < decoded[decodedIndex].size())
                {
                    k++;
                    it++;
                }
                
                //fingerprint already exists, so add to it
                if(compareRemainders(r, it->remainder, this->r))
                {
                    decoded[decodedIndex][k].count++;
                }
                else if(k + 1 == decoded[decodedIndex].size())
                {
                    //cout << "two" << endl;
                    memcpy(temp.remainder, r, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
                    //temp.remainder = r;
                    temp.count = 1;
                    if(r[0] < it->remainder[0])
                    {
                        decoded[decodedIndex].insert(it, temp);
                    }
                    else
                    {
                        decoded[decodedIndex].push_back(temp);
                    }
                }
                else//fingerprint does not exist so insert fingerprint in correct location to maintain ascending order
                {
                    //cout << "three" << endl;
                    memcpy(temp.remainder, r, sizeof(unsigned char) * (this->r / BITS_PER_UCHAR));
                    //temp.remainder = r;
                    temp.count = 1;
                    decoded[decodedIndex].insert(it, temp);
                }
            }
            if(toEncode)
            {
                encode(info);
            }
            clearDecodedStructure(info.decodeSize);
        }
    }
    
    
    
    unsigned long getLongFromChar(unsigned char binaryKmer[K_SIZE])
    {
        unsigned long toRet = 0x00000000;
        
        for(int i = 0; i < K_SIZE; i++)
        {
            toRet |= (long)binaryKmer[i] << (((K_SIZE - 1) - i) * 8);
        }
        
        return toRet;
    }
    
    unsigned long getLongFromChar(unsigned char *binaryKmer, int sizeR)
    {
        unsigned long toRet = 0x00000000;
        int numIteration = sizeR / BITS_PER_UCHAR;
        
        for(int i = 0; i < numIteration; i++)
        {
            toRet |= (long)binaryKmer[i] << (((numIteration - 1) - i) * 8);
        }
        
        return toRet;
    }
    
    
    
    
    
    
public:
    quotientFilter(int sizeQ, int sizeR)
    {
        long m = pow(2, sizeQ);
        
        this->r = sizeR;
        this->q = sizeQ;
        this->size = m + 100;
        
        this->tab_entry = new struct QF_Entry[m + 100];
        
        cout << "QF has: " << m << " entries" << endl;
        
        for(int i = 0; i < m; i++)
        {
            this->tab_entry[i].bitHolder = 0x00;
            this->tab_entry[i].count = 0;
        }
        
        this->initDecodedStructure();
    }
    
    void setKmerUsed(char kmer[WIN_SIZE_PLUS])
    {
        unsigned char result[K_SIZE], *remainder = NULL;
        unsigned long hash;
        memset(result, 0x00, sizeof(unsigned char)*K_SIZE);
        this->convert_char_to_binary(kmer, result);
        hash = this->convert_bin_to_hash(result, this->q);
        remainder = this->convert_bin_to_remainder(result, this->q, this->r);
        int distCount = -1;
        
        //if this location is shifted, it means we have to decode it
        if(this->tab_entry[hash].bitHolder & IS_SHIFTED)
        {
            //cout << "this" << endl;
            struct decodeInfo info = this->decodeCluster(hash);
            unsigned long decodedIndex = hash - info.clusterStartLocation;
            
            
            for(int i = 0; i < decodedIndex; i++)
            {
                distCount += decoded[i].size();
            }
            for(int i = 0; i < decoded[decodedIndex].size(); i++)
            {
                distCount++;
                if(this->compareRemainders(remainder, decoded[decodedIndex][i].remainder, this->r))
                {
                    //clearDecodedStructure(info.decodeSize);
                    this->tab_entry[info.clusterStartLocation + distCount].bitHolder |= USED;
                    break;
                }
            }
            
            //make sure to clear the decoded structure
            clearDecodedStructure(info.decodeSize);
        }
        else
        {
            if(this->compareRemainders(remainder, this->tab_entry[hash].remainder, this->r))
            {
                this->tab_entry[hash].bitHolder |= USED;
            }
            else
            {
                hash++;
                while((this->tab_entry[hash].bitHolder & IS_CONTINUATION) && hash < this->size)
                {
                    if(this->compareRemainders(remainder, this->tab_entry[hash].remainder, this->r))
                    {
                        this->tab_entry[hash].bitHolder |= USED;
                    }
                    hash++;
                }
            }
        }
    }
    
    short getKmerCount(char kmer[WIN_SIZE_PLUS], bool &used)
    {
        unsigned char result[K_SIZE], *remainder = NULL;
        unsigned long hash;
        memset(result, 0x00, sizeof(unsigned char)*K_SIZE);
        this->convert_char_to_binary(kmer, result);
        hash = this->convert_bin_to_hash(result, this->q);
        remainder = this->convert_bin_to_remainder(result, this->q, this->r);
        short count;
        
        //if this location is shifted, it means we have to decode it
        if(this->tab_entry[hash].bitHolder & IS_SHIFTED)
        {
            //cout << "this" << endl;
            struct decodeInfo info = this->decodeCluster(hash);
            unsigned long decodedIndex = hash - info.clusterStartLocation;
            
            for(int i = 0; i < decoded[decodedIndex].size(); i++)
            {
                if(this->compareRemainders(remainder, decoded[decodedIndex][i].remainder, this->r))
                {
                    if(decoded[decodedIndex][i].used)
                    {
                        used = true;
                    }
                    else
                    {
                        used = false;
                    }
                    count = decoded[decodedIndex][i].count;
                    clearDecodedStructure(info.decodeSize);
                    return count;
                }
            }
            
            //make sure to clear the decoded structure
            clearDecodedStructure(info.decodeSize);
        }
        else
        {
            if(this->tab_entry[hash].count == 0)
            {
                //do nothing
            }
            else if(this->compareRemainders(remainder, this->tab_entry[hash].remainder, this->r))
            {
                if(this->tab_entry[hash].bitHolder & USED)
                {
                    used = true;
                }
                else
                {
                    used = false;
                }
                return this->tab_entry[hash].count;
            }
            else
            {
                hash++;
                while((this->tab_entry[hash].bitHolder & IS_CONTINUATION) && hash < this->size)
                {
                    if(this->compareRemainders(remainder, this->tab_entry[hash].remainder, this->r))
                    {
                        if(this->tab_entry[hash].bitHolder & USED)
                        {
                            used = true;
                        }
                        else
                        {
                            used = false;
                        }
                        return this->tab_entry[hash].count;
                    }
                    hash++;
                }
            }
        }
        //this means none exist
        used = false;
        return 0;
    }
    
    void insertKmer(char kmer[WIN_SIZE_PLUS])
    {
        unsigned char result[K_SIZE], *remainder;
        unsigned long hash;
        
        memset(result, 0x00, sizeof(unsigned char)*K_SIZE);
        this->convert_char_to_binary(kmer, result);
        hash = this->convert_bin_to_hash(result, this->q);
        remainder = this->convert_bin_to_remainder(result, this->q, this->r);
        
        this->insertElement(hash, remainder, result);
    }
    
    void convert_char_to_binary (char *kmer_name, unsigned char *result)
    {
        
        int i=0,j=0,k=0;
        
        for (i=0; i<WIN_SIZE; i++)
        {
            if ((i%4 == 0) && (i>0))
            {
                j++;
                k=0;
            }
            
            switch (kmer_name[i])
            {
                case 'A':
                    result[j] |= 0 << (7 - (2*k));
                    result[j] |= 0 << (6 - (2*k));
                    break;
                    
                case 'C':
                    result[j] |= 0 << (7 - (2*k));
                    result[j] |= 1 << (6 - (2*k));
                    break;
                    
                case 'G':
                    result[j] |= 1 << (7 - (2*k));
                    result[j] |= 0 << (6 - (2*k));
                    break;
                    
                case 'T':
                    result[j] |= 1 << (7 - (2*k));
                    result[j] |= 1 << (6 - (2*k));
                    break;
            }
            k++;
        }
    }
    
    void convert_binary_to_char (unsigned char *result, char output_str[WIN_SIZE_PLUS])
    {
        
        int i=0,j=0,len=0;
        //len = *out_len;
        
        for (i=0; i<K_SIZE; i++)
        {
            for (j=7; j>=0; j-=2)
            {
                unsigned mask1 = createMask(j-1,j);
                unsigned temp = mask1 & result[i];
                if (j>1)
                    temp = temp >> (j-1);
                //printf("i: %d, j: %d, mask1: %d, input: %d, temp: %d \n", i, j, mask1, result[i], temp);
                
                switch (temp)
                {
                    case 0:
                        strcpy(&output_str[len],"A");
                        break;
                        
                    case 1:
                        strcpy(&output_str[len],"C");
                        break;
                        
                    case 2:
                        strcpy(&output_str[len],"G");
                        break;
                        
                    case 3:
                        strcpy(&output_str[len],"T");
                        break;
                }
                len++;
            }     
        }
        
        output_str[len] = '\0';
        //*out_len = len;
        
    }
    
    string convert_binary_to_string(long num, int numBinaryDigits)
    {
        string toRet = "";
        long temp;
        
        for(int i = 0; i < numBinaryDigits; i += 2)
        {
            temp = num & 3;
            switch (temp) {
                case 0:
                    toRet.insert(0,"A");
                    break;
                    
                case 1:
                    toRet.insert(0,"C");
                    break;
                    
                case 2:
                    toRet.insert(0,"G");
                    break;
                    
                case 3:
                    toRet.insert(0,"T");
                    break;
            }
            
            num = num >> 2;
        }
        
        return toRet;
    }
    
    string runQR(long q, unsigned char *r)
    {
        int s_pointer = (this->q / 2), lastChar = WIN_SIZE - 1;
        int S_Size_binary = (WIN_SIZE * 2) - (this->q + this->r);
        int S_Size_chars = S_Size_binary / 2, s_loc = s_pointer - 1;
        string contig = "", temp, maxKmer;
        char testChars[] = "ACGT", kmer[WIN_SIZE_PLUS] = "";
        unsigned short maxFreq = 0, freq = 0;
        bool isKmerUsed;
        //initiate the contig with the q
        contig = this->convert_binary_to_string(q, this->q);
        //fill our S with all 'A's
        for(int j = 0; j < S_Size_chars; j++)
        {
            contig += "A";
        }
        //printf("r0: %d, r1: %d\n", r[0], r[1]);
        //fill the remainder of our start contig with the r
        unsigned long longR = this->getLongFromChar(r, this->r);
        contig += this->convert_binary_to_string(longR, this->r);

        //cout << "Start: " << contig << endl;
        //temp = contig;
        //if(this->is_start(contig))
        //{
            for(int j = 0; j < S_Size_chars; j++)
            {
                temp = contig.substr(j + 1);
                for(int k = 0; k < 4; k++)
                {
                    temp[s_loc] = testChars[k];
                    for(int n = 0; n < 4; n++)
                    {
                        if(temp.size() < 32)
                        {
                            temp += testChars[n];
                        }
                        else
                        {
                            temp[lastChar] = testChars[n];
                        }
                    
                        //cout << "trying: " << temp << endl;
                        strcpy(kmer, temp.c_str());
                        freq = this->getKmerCount(kmer, isKmerUsed);
                    
                        if(!isKmerUsed)
                        {
                            if(freq > maxFreq)
                            {
                                maxFreq = freq;
                                maxKmer = temp;
                            }
                        }
                    }
                }
                
                if(maxFreq == 0)
                {
                    contig = "";
                    return contig;
                }
            
                //cout << "Found: " << maxKmer << endl;
                //cout << "Freq:  " << maxFreq << endl;
                contig[s_pointer] = maxKmer[s_loc];
                contig += maxKmer[lastChar];
                strcpy(kmer, maxKmer.c_str());
                setKmerUsed(kmer);
                //cout << "cur:   " << contig << endl;
                s_pointer++;
                maxFreq = 0;
                maxKmer = "";
            }
        
        maxFreq = 1;
        while(maxFreq != 0)
        {
            maxFreq = 0;
            maxKmer = "";
            //temp should now be 31 characters long
            temp = contig.substr(contig.size() - WIN_SIZE + 1);
            
            for(int k = 0; k < 4; k++)
            {
                if(temp.size() < WIN_SIZE)
                {
                    temp += testChars[k];
                }
                else
                {
                    temp[lastChar] = testChars[k];
                }
                
                strcpy(kmer, temp.c_str());
                freq = this->getKmerCount(kmer, isKmerUsed);
                
                if(!isKmerUsed)
                {
                    if(freq > maxFreq)
                    {
                        maxFreq = freq;
                        maxKmer = temp;
                    }
                }
                
            }
            
            if(maxFreq != 0)
            {
                strcpy(kmer, maxKmer.c_str());
                setKmerUsed(kmer);
                contig += maxKmer[lastChar];
            }
        }
        
            return contig;
        //}
        //else
        //{
        //    contig = "";
        //    return contig;
        //}
    }
    
    bool is_start(string startContig)
    {
        char testChar[5] = "ACGT", kmer[WIN_SIZE_PLUS] = "";;
        //cout << "Before: " << startContig << endl;
        startContig.pop_back();
        startContig.insert(0, "A");
        //cout << "After:  " << startContig << endl;
        for(int i = 0; i < 4; i++)
        {
            startContig[0] = testChar[i];
            for(int k = 0; k < 4; k++)
            {
                startContig[28] = testChar[k];
                strcpy(kmer, startContig.c_str());
                //if(this->getKmerCount(kmer))
                //{
                    return false;
                //}
            }
        }
        
        return true;
    }
    
    void build_S()
    {
        long i = 0,f, tempHash;
        string contig = "", temp, maxKmer, returnContig;
        unsigned char tempRemainder[2];
        struct decodeInfo info;
        vector<string> contigList;
        ofstream outfile("myContig.txt");
        int decodeCounter = 0;
        
        while(i < this->size)//this->tab_entry[i].count < 300 || (this->tab_entry[i].bitHolder & IS_SHIFTED))
        {
            if(this->tab_entry[i].bitHolder & IS_SHIFTED)
            {
                info = this->decodeCluster(i);
                decodeCounter = 0;
                for(f = 0; f < info.decodeSize; f++)
                {
                    for(int l = 0; l < decoded[f].size(); l++)
                    {
                        if(decoded[f][l].count >= 20)
                        {
                            if(!(this->tab_entry[i + decodeCounter - 1].bitHolder & USED))
                            {
                                //set the new starting QR as used
                                this->tab_entry[i + decodeCounter - 1].bitHolder |= USED;
                                //printf("tab: %d, decode: %d\n", this->tab_entry[i+ decodeCounter - 1].remainder[0], decoded[f][l].remainder[0]);
                                memcpy(tempRemainder, decoded[f][l].remainder, sizeof(unsigned char) * (this->r / 8));
                                this->clearDecodedStructure(info.decodeSize);
                                //we subtract one because i is actually pointing to the spot just after the start of the cluster
                                tempHash = i+f-1;
                                returnContig = runQR(tempHash, tempRemainder);
                                if(returnContig.length() > 1)
                                {
                                    contigList.push_back(returnContig);
                                }
                                info = this->decodeCluster(i);
                            }
                        }
                        decodeCounter++;
                    }
                }
                i += info.clusterSize - 1;
                this->clearDecodedStructure(info.decodeSize);
            }
            else if(this->tab_entry[i].count > 40)
            {
                if(!(this->tab_entry[i].bitHolder & USED))
                {
                    this->tab_entry[i].bitHolder |= USED;
                    returnContig = runQR(i, this->tab_entry[i].remainder);
                    if(returnContig.length() > 1)
                    {
                        contigList.push_back(returnContig);
                    }
                }
            }
            
            
            i++;
        }
        
        int contigCounter = 1;
        for(auto thing : contigList)
        {
            outfile << ">contig_" << contigCounter << endl;
            outfile << thing << endl;
            contigCounter++;
        }
        outfile.close();
        /*
        cout << "count: " << this->tab_entry[i].count << endl;
        
        contig = this->convert_binary_to_string(i, this->q);
        for(int j = 0; j < S_Size_chars; j++)
        {
            contig += "A";
        }
        //printf("r: %d\n", this->tab_entry[i].remainder[0]);
        contig += this->convert_binary_to_string((long)this->tab_entry[i].remainder[0], this->r);
        //1101101001110100
        cout << "Start: " << contig << endl;
        
        for(int j = 0; j < S_Size_chars; j++)
        {
            temp = contig.substr(j + 1);
            for(int k = 0; k < 4; k++)
            {
                temp[s_loc] = testChars[k];
                for(int n = 0; n < 4; n++)
                {
                    if(temp.size() < 32)
                    {
                        temp += testChars[n];
                    }
                    else
                    {
                        temp[lastChar] = testChars[n];
                    }
                    
                    //temp = "GACCTGCTGTACCGCAACTCGTGGAACGAAGT";
                    cout << "trying: " << temp << endl;
                    strcpy(kmer, temp.c_str());
                    freq = this->getKmerCount(kmer);
                    
                    if(freq > maxFreq)
                    {
                        maxFreq = freq;
                        maxKmer = temp;
                    }
                }
            }
            
            cout << "Found: " << maxKmer << endl;
            cout << "Freq:  " << maxFreq << endl;
            contig[s_pointer] = maxKmer[s_loc];
            contig += maxKmer[lastChar];
            cout << "cur:   " << contig << endl;
            s_pointer++;
            maxFreq = 0;
            maxKmer = "";
        }
        
        
        
        cout << "built: " << contig << endl;*/
        
        
        
    }
    
};




#endif
