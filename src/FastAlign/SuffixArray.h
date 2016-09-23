#ifndef _SUFFIX_ARRAY_H
#define _SUFFIX_ARRAY_H

#if defined(OPEN_MP)
  #include <parallel/algorithm>
#else
  #include <algorithm>
#endif


#include <map>
#include <string>
#include <sstream>
#include <stdint.h>
#include "base/SVector.h"
#include "src/DNAVector.h"
#include "base/FileParser.h"
#include "extern/logger/log.h"
#include "src/Cola/Cola.h"
#include "src/FastAlign/AlignmentParams.h"
#include "src/FastAlign/SeedingObjects.h"
#include "src/FastAlign/DNASeqs.h"


//======================================================
/** Suffix Array Element */
class SuffixArrayElement {
public:
    SuffixArrayElement(): m_index(), m_offset(), m_strand(1) {}
    SuffixArrayElement(int i, int o, int s): m_index(i), m_offset(o), m_strand(s) {}

    bool operator<(const SuffixArrayElement& other) const { 
        if(m_index!=other.m_index) {
            return (m_index<other.m_index);
        } else {
            return (m_offset<other.m_offset);
        }  
    }

    int getIndex() const    { return m_index;  }
    int getOffset() const   { return m_offset; }
    int getStrand() const   { return m_strand; }

    string toString() const;

private:
    int  m_index;   /// Index of string to which the substring comes from
    int  m_offset;  /// Offset into the string which the substring has been extracted from
    int  m_strand;  /// true: positive strand  false:negative strand   //TODO this should only be a bool but int to adapt to Search code
};

//======================================================

//======================================================
/** Object containing actual sequence reference and offset 
    This is to be used for comparison purpose when sequence data is
    not the same as suffix data */
template<class StringType>
class SequenceWithOffset {
public:
    SequenceWithOffset(const StringType& sequence, int offset): m_sequence(sequence), m_offset(offset) {}

    const StringType&  getSequence() const { return m_sequence; } 
    int getOffset() const                  { return m_offset;   }

private:
    const StringType& m_sequence;   /// Sequence content 
    int  m_offset;            /// Offset into the string which the substring has been extracted from
};
//======================================================

//======================================================
template<class StringContainerType, class StringType>
class SuffixArray {
public:
    // Ctor1:
    SuffixArray(const StringContainerType& strings, int stepSize): m_suffixes(), m_strings(strings),
                m_stepSize_p(stepSize) { 
      constructSuffixes(); 
    }

    const StringContainerType& getSeqs() { return m_strings; }

    int   getSuffixStep() const                             { return m_stepSize_p;                 }  

    int getSize() const                                     { return m_suffixes.size();            } 
    SuffixArrayElement getByIndex(unsigned long idx) const  { return m_suffixes[idx];              } 
    const svec<SuffixArrayElement>& getSuffixes() const     { return m_suffixes;                   }

    void constructSuffixes(); 
    void sortSuffixes(); 
    bool filterLowComp(const SuffixArrayElement& sr);  //Determine if a substring is of very low complexity (repeat character)

    int getStringSize(int sIdx) const                       { return m_strings[sIdx].isize();      }
    const StringType& getString(int sIdx) const             { return m_strings[sIdx];              }  

    string getSeq(const SuffixArrayElement& sr, int startIdx, int endIdx) const; 
    string getSeq(const SuffixArrayElement& sr, int len) const; 
    string getSeq(const SuffixArrayElement& sr) const; 
    string getSeq(unsigned long index) const;  //Get Substring sequence chars via the index in substrings
 
    void getDNA(const SuffixArrayElement& sr, int startIdx, int endIdx, StringType& outDNA)  const; 

    int compareBases(const SuffixArrayElement& s1, const SuffixArrayElement& s2) const;
    int compareBases(const SuffixArrayElement& s1, const StringType& d2) const;
    int compareBases(const SuffixArrayElement& s1, const StringType& d2, int offset2) const; 

    struct CmpSuffixArrayElementOL { //SuffixArrayElement comparison struct for overlap finding
        CmpSuffixArrayElementOL(const SuffixArray<StringContainerType, StringType>& s) : m_substrings(s) {}
        bool operator() (const SuffixArrayElement& s1, const SuffixArrayElement& s2) const { 
            return (m_substrings.compareBases(s1, s2) == -1);
        }
        bool operator() (const SuffixArrayElement& s1, const StringType& d2) const { 
            return (m_substrings.compareBases(s1, d2) == -1);
        }
        bool operator() (const SuffixArrayElement& s1, const SequenceWithOffset<StringType>& swo) const { 
            return (m_substrings.compareBases(s1, swo.getSequence(), swo.getOffset()) == -1);
        }
        const SuffixArray<StringContainerType, StringType>& m_substrings;
    };


private:
    void  getSeq(const SuffixArrayElement& sr, string& outSeq) const    { outSeq = getSeq(sr); }

    svec<SuffixArrayElement>      m_suffixes;     /// Vector of suffixes 
    const StringContainerType&    m_strings;      /// Reference to the list of strings from which suffixes where constructed
    int                           m_stepSize_p;   /// Parameter specifing the size of steps for constructing suffixes 
};
//======================================================

//======================================================
template<class StringContainerType, class StringType>
void SuffixArray<StringContainerType, StringType>::getDNA(const SuffixArrayElement& sr, int startIdx, int endIdx, StringType& outDNA) const {  
    int stringSize = m_strings[sr.getIndex()].size();
    if(endIdx>=stringSize) { endIdx = stringSize-1; }
    int len = endIdx-startIdx+1;
    if(sr.getStrand()==1) {
        outDNA.SetToSubOf(m_strings.getSeqByIndex(sr.getIndex()), startIdx, len); 
    } else {
        outDNA.SetToSubOf(m_strings.getSeqRCByIndex(sr.getIndex()), startIdx, len);
    }

}  

template<class StringContainerType, class StringType>
string SuffixArray<StringContainerType, StringType>::getSeq(const SuffixArrayElement& sr, int startIdx, int endIdx) const { 
    int stringSize = m_strings[sr.getIndex()].size();
    if(endIdx>=stringSize) { endIdx = stringSize-1; }
    int to = endIdx-startIdx+1;
    if(sr.getStrand()==1) {
        return m_strings.getSeqByIndex(sr.getIndex()).Substring(startIdx, to); 
    } else {
        return m_strings.getSeqRCByIndex(sr.getIndex()).Substring(startIdx, to); 
    }
}  

template<class StringContainerType, class StringType>
string SuffixArray<StringContainerType, StringType>::getSeq(const SuffixArrayElement& sr, int len) const { 
    //TODO warn if given length exceeds the available sequence
    return getSeq(sr, sr.getOffset(), sr.getOffset()+len-1);
}

template<class StringContainerType, class StringType>
string SuffixArray<StringContainerType, StringType>::getSeq(const SuffixArrayElement& sr) const { 
    return getSeq(sr, sr.getOffset(), m_strings[sr.getIndex()].size()-1);
}

template<class StringContainerType, class StringType>
string SuffixArray<StringContainerType, StringType>::getSeq(unsigned long index) const { 
    SuffixArrayElement sr = m_suffixes[index];
    return getSeq(sr, sr.getOffset(), m_strings[sr.getIndex()].size()-1);
}

template<class StringContainerType, class StringType>
void SuffixArray<StringContainerType, StringType>::constructSuffixes() {
    m_suffixes.clear();
    if(m_strings.getNumSeqs()==0) { return; } //There are no input strings to continue with
    double numSubstrings = (double) m_strings.getNumSeqs()*m_strings.getSize(0)/getSuffixStep()+1;
    m_suffixes.reserve(numSubstrings); // If strings have varying sizes, this is just an estimate
    for(int i=0; i<m_strings.getNumSeqs(); i++) {
        for(int j=0; j<m_strings[i].isize(); j+=getSuffixStep()) {
            SuffixArrayElement sr(i, j, 1); //i:index j:offset 
            if(filterLowComp(sr)) { continue; }
            m_suffixes.push_back(sr); //i:index j:offset 
        }
    }
    sortSuffixes();
    FILE_LOG(logDEBUG) <<"Total number of strings: " << m_strings.getNumSeqs();
    cout <<"Total number of strings: " << m_strings.getNumSeqs() << endl;
    FILE_LOG(logDEBUG) <<"Total number of substrings: " << m_suffixes.size();
    cout <<"Total number of substrings: " << m_suffixes.size() << endl;
} 

template<class StringContainerType, class StringType>
bool SuffixArray<StringContainerType, StringType>::filterLowComp(const SuffixArrayElement& sr) {
    int  idx    = sr.getIndex();
    int  offset = sr.getOffset();
    int  end    = min(offset + getSuffixStep(), m_strings[idx].isize());
    char prevC  = ' '; // Previous rep character
    char currC  = ' '; // Current character
    int  repCnt = 0;   // keep track of size of longest repeat contig
    int  nCnt   = 0;   // Keep track of nonesense characters
    for(int i=offset; i<end; i++) {
        char currC =  m_strings[idx][i];
        if(currC=='n' || currC=='N') { nCnt++; }
        if(prevC==currC) {
            repCnt++;
        } else {
            repCnt = 0;
            prevC  = currC;
        }
    }
    if((float)repCnt/getSuffixStep()>0.90 || nCnt>0.5*getSuffixStep()) { 
        return true;
    } else {
        return false;
    }
}

template<class StringContainerType, class StringType>
void SuffixArray<StringContainerType, StringType>::sortSuffixes() {
    FILE_LOG(logDEBUG) << "Starting to sort SuffixArray";
    cout << "Starting to sort SuffixArray" << endl;
    // Could not do this in header by including namespaces as std has been declared in a parent header.
    #if defined(OPEN_MP)
        __gnu_parallel::stable_sort(m_suffixes.begin(), m_suffixes.end(), CmpSuffixArrayElementOL(*this));
    #else
        std::stable_sort(m_suffixes.begin(), m_suffixes.end(), CmpSuffixArrayElementOL(*this));
    #endif
    FILE_LOG(logDEBUG) << "Finished sorting SuffixArray";
    cout << "Finished sorting suffixes" << endl;
}

template<class StringContainerType, class StringType>
int SuffixArray<StringContainerType, StringType>::compareBases(const SuffixArrayElement& s1, const SuffixArrayElement& s2) const { 
    int idx1    = s1.getIndex();
    int offset1 = s1.getOffset();
    int strand1 = s1.getStrand();
    int idx2    = s2.getIndex();
    int offset2 = s2.getOffset();
    int strand2 = s2.getStrand();

    int limit = min(getStringSize(idx1)-offset1, getStringSize(idx2)-offset2);

    const StringType& d1 = getString(idx1);
    const StringType& d2 = getString(idx2); 

    for(int i=0; i<limit; i++) {
        if(d1[offset1+i]<d2[offset2+i]) {
            return -1; //Smaller
        } else if(d1[offset1+i]>d2[offset2+i]) {
            return 1;  //Larger
        }
    }
    return 0; // Same
}

template<class StringContainerType, class StringType>
int SuffixArray<StringContainerType, StringType>::compareBases(const SuffixArrayElement& s1, const StringType& d2) const { 
    return compareBases(s1, d2, 0); //Offset of d2 is 0, as the actual sequence is passed here
}

template<class StringContainerType, class StringType>
int SuffixArray<StringContainerType, StringType>::compareBases(const SuffixArrayElement& s1, const StringType& d2, int offset2) const { 
    int idx1    = s1.getIndex();
    int offset1 = s1.getOffset();
    int strand1 = s1.getStrand();

    int limit = min(getStringSize(idx1)-offset1, d2.size()-offset2);
    const StringType& d1 = getString(idx1);
    for(int i=0; i<limit; i++) {
        if(d1[offset1+i]<d2[offset2+i]) {
            return -1; //Smaller
        } else if(d1[offset1+i]>d2[offset2+i]) {
            return 1;  //Larger
        }
    }
    return 0; // Same
}
 
 //======================================================
#endif //_SUFFIX_ARRAY_H
