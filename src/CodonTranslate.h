
#ifndef CODONTRANSLATE_H
#define CODONTRANSLATE_H

#include "base/SVector.h"
#include "src/DNAVector.h"
#include <string>
#include <map>

inline char as_i(char b) {
  char i = -1;

  switch(b) {
  case 'A':
    i = 0;
    break;
  case 'C':
    i = 1;
    break;
  case 'G':
    i = 2;
    break;
  case 'T':
    i = 3;
    break;
  default:   
    break;
  }
  return i;
}


class CodonMatrix
{
 public:
  CodonMatrix() {
    m_row.resize(21, 0);
  }
  void Set(int i, double v) {
    m_row[i] = v;
  }

  double Get(int i) const {return m_row[i];}
  int Size() const {return m_row.isize();}
 private:
  svec<double> m_row;
};


class CodonTranslator
{
public:
  CodonTranslator() {
    init();
  }

  void init(); 

/*
  int GetCodonIndexByIndex(char one, char two, char three) {
    int i = GetIndex(one, two, three);
    return m_index[i];
  }
  
  char GetCodonByIndex(char one, char two, char three) {
    int i = GetIndex(one, two, three);
    return m_codons[m_index[i]];
  }
*/
  int GetCodonIndex(char one, char two, char three) const {
    int i = GetIndex(as_i(one), as_i(two), as_i(three));
    return m_index[i];
  }
  
  char GetCodon(char one, char two, char three) const {
    int i = GetIndex(as_i(one), as_i(two), as_i(three));
    if(i>m_index.isize()) { return 'N'; }
    return m_codons[m_index[i]];
  }

  const CodonMatrix & GetMatrix(int idx) const {return m_matrix[idx];}
 
  /** Outputs synonymous mutation analysis for a give sequence */
  void GetSynMutInfo(const DNAVector& seq, string& outStr); 

private:
  /** Returns codons corresponding to a given aminoacid character */
  const svec<string>& GetCodons(char aa) const { 
    map<char, svec<string> >::const_iterator pos = m_amino2codons.find(aa);
    if(pos==m_amino2codons.end()) {
      pos = m_amino2codons.find('N'); // Place holder that has been created on initiation of class
    }
    return pos->second; 
  }
  /** Returns codons corresponding and equivalet to a given codon, including the input codon itself */
  const svec<string>& GetEqCodons(const string& aa) const { return GetCodons(GetCodon(aa[0], aa[1], aa[2])); }
  int GetIndex(char one, char two, char three) const;
  void Set(const char * codon, const char * three, int index); 
  void SetFullMatrix();
  void SetMatrix(char a, char b, double v); 
  int Find(char a) const; 

  svec<char> m_codons;
  svec<int> m_index;                       /// Index of codons corresponding to integer value of codon 
  map<char, svec<string> > m_amino2codons; /// Mapping aminoacid characters to the equivalent codons
  svec<CodonMatrix> m_matrix;
  map<string, string> eqCodonLTable;       /// Equivalent codon lookup table (filled only when needed)
};

#endif //CODONTRANSLATE_H

