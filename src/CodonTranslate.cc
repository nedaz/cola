#include "src/CodonTranslate.h"
#include "extern/logger/log.h"

int CodonTranslator::GetIndex(char one, char two, char three) const {
  
  int i = one;
  i *= 4;
  i += (int)two;
  i *= 4;
  i += (int)three;
  return i;
}

void CodonTranslator::Set(const char * codon, const char * three, int index) {
  
  int i = GetIndex(as_i(three[0]), as_i(three[1]), as_i(three[2]));
  m_index[i] = index;
  //cout << "i=" << i << " index=" << index << endl;
  m_codons[index] = codon[0];
  m_amino2codons[codon[0]].push_back(three);
}

void CodonTranslator::SetMatrix(char a, char b, double v) {
  int i;
  int idx1 = Find(a);
  int idx2 = Find(b);
  if (idx1 == -1 || idx2 == -1)
    return;

  CodonMatrix & m = m_matrix[idx1];
  m.Set(idx2, v);
}

int CodonTranslator::Find(char a) const {
  for (int i=0; i<m_codons.isize(); i++) {
    if (m_codons[i] == a) {
      return i;
    }
  }
  return -1;
}

void CodonTranslator::GetSynMutInfo(const DNAVector& seq, string& outStr) {
  // Create output lookup table
  for(int i=0; i<seq.size()-2; i++) { 
    string codon = seq.Substring(i, 3);
    FILE_LOG(logDEBUG4) << codon; 
    string lookup = eqCodonLTable[codon];
    if(lookup == "") {
      svec<string> eqCodons = GetEqCodons(codon);
      string pos0, pos1, pos2;
      for(int c=0; c<eqCodons.isize(); c++) {
        string eqCodon = eqCodons[c];
        FILE_LOG(logDEBUG4) << eqCodon; 
        if(eqCodon==codon) { continue; }
        if(eqCodon[1]==codon[1] && eqCodon[2]==codon[2]) {
          pos0 += eqCodon[0];
        }
        if(eqCodon[0]==codon[0] && eqCodon[2]==codon[2]) {
          pos1 += eqCodon[1];
        }
        if(eqCodon[0]==codon[0] && eqCodon[1]==codon[1]) {
          pos2 += eqCodon[2];
        }
      }
      eqCodonLTable[codon] = pos0 + "," + pos1 + "," + pos2;
    }
    lookup = eqCodonLTable[codon];
    FILE_LOG(logDEBUG4) << lookup; 
    outStr += (lookup+","); 
  }
  
}
 
void CodonTranslator::init() {
  m_index.resize(64, 0);
  m_codons.resize(21, 0);
  
   m_amino2codons['N'].push_back(""); //Unknown Codon

  Set("A", "GCT", 0); 
  Set("A", "GCC", 0); 
  Set("A", "GCA", 0); 
  Set("A", "GCG", 0); 

  Set("L", "TTA", 1); 
  Set("L", "TTG", 1); 
  Set("L", "CTT", 1); 
  Set("L", "CTC", 1); 
  Set("L", "CTA", 1); 
  Set("L", "CTG", 1);

  Set("R", "CGT", 2); 
  Set("R", "CGC", 2); 
  Set("R", "CGA", 2); 
  Set("R", "CGG", 2); 
  Set("R", "AGA", 2); 
  Set("R", "AGG", 2);
 
  Set("K", "AAA", 3); 
  Set("K", "AAG", 3); 

  Set("N", "AAT", 4); 
  Set("N", "AAC", 4); 

  Set("M", "ATG", 5); 

  Set("D", "GAT", 6); 
  Set("D", "GAC", 6); 

  Set("F", "TTT", 7); 
  Set("F", "TTC", 7); 

  Set("C", "TGT", 8); 
  Set("C", "TGC", 8); 

  Set("P", "CCT", 9); 
  Set("P", "CCC", 9); 
  Set("P", "CCA", 9); 
  Set("P", "CCG", 9); 

  Set("Q", "CAA", 10); 
  Set("Q", "CAG", 10); 

  Set("S", "TCT", 11); 
  Set("S", "TCC", 11); 
  Set("S", "TCA", 11); 
  Set("S", "TCG", 11); 
  Set("S", "AGT", 11); 
  Set("S", "AGC", 11); 

  Set("E", "GAA", 12); 
  Set("E", "GAG", 12); 

  Set("T", "ACT", 13); 
  Set("T", "ACC", 13); 
  Set("T", "ACA", 13); 
  Set("T", "ACG", 13); 

  Set("G", "GGT", 14); 
  Set("G", "GGC", 14); 
  Set("G", "GGA", 14); 
  Set("G", "GGG", 14); 

  Set("W", "TGG", 15); 

  Set("H", "CAT", 16); 
  Set("H", "CAC", 16); 

  Set("Y", "TAT", 17); 
  Set("Y", "TAC", 17); 


  Set("I", "ATT", 18); 
  Set("I", "ATC", 18); 
  Set("I", "ATA", 18); 

  Set("V", "GTT", 19); 
  Set("V", "GTC", 19); 
  Set("V", "GTA", 19); 
  Set("V", "GTG", 19); 

  Set("*", "TAG", 20); 
  Set("*", "TGA", 20); 
  Set("*", "TAA", 20); 

  m_matrix.resize(21);
    
  SetFullMatrix();
}

void CodonTranslator::SetFullMatrix() {
  SetMatrix('A', 'A', 0.21558);
  SetMatrix('A', 'R', -0.00181159);
  SetMatrix('A', 'N', -0.0452899);
  SetMatrix('A', 'D', -0.0452899);
  SetMatrix('A', 'C', 0.0416667);
  SetMatrix('A', 'Q', -0.00181159);
  SetMatrix('A', 'E', -0.00181159);
  SetMatrix('A', 'G', 0.0416667);
  SetMatrix('A', 'H', -0.0452899);
  SetMatrix('A', 'I', -0.00181159);
  SetMatrix('A', 'L', -0.00181159);
  SetMatrix('A', 'K', -0.00181159);
  SetMatrix('A', 'M', -0.00181159);
  SetMatrix('A', 'F', -0.0452899);
  SetMatrix('A', 'P', -0.00181159);
  SetMatrix('A', 'S', 0.0851449);
  SetMatrix('A', 'T', 0.0416667);
  SetMatrix('A', 'W', -0.0887681);
  SetMatrix('A', 'Y', -0.0452899);
  SetMatrix('A', 'V', 0.0416667);
  SetMatrix('A', 'B', -0.0452899);
  SetMatrix('A', 'Z', -0.00181159);
  SetMatrix('A', 'X', 0.0416667);
  SetMatrix('A', '*', -0.132246);
  SetMatrix('R', 'A', 0.00462963);
  SetMatrix('R', 'R', 0.226852);
  SetMatrix('R', 'N', 0.0416667);
  SetMatrix('R', 'D', -0.0324074);
  SetMatrix('R', 'C', -0.0694444);
  SetMatrix('R', 'Q', 0.0787037);
  SetMatrix('R', 'E', 0.0416667);
  SetMatrix('R', 'G', -0.0324074);
  SetMatrix('R', 'H', 0.0416667);
  SetMatrix('R', 'I', -0.0694444);
  SetMatrix('R', 'L', -0.0324074);
  SetMatrix('R', 'K', 0.115741);
  SetMatrix('R', 'M', 0.00462963);
  SetMatrix('R', 'F', -0.0694444);
  SetMatrix('R', 'P', -0.0324074);
  SetMatrix('R', 'S', 0.00462963);
  SetMatrix('R', 'T', 0.00462963);
  SetMatrix('R', 'W', -0.0694444);
  SetMatrix('R', 'Y', -0.0324074);
  SetMatrix('R', 'V', -0.0694444);
  SetMatrix('R', 'B', 0.00462963);
  SetMatrix('R', 'Z', 0.0416667);
  SetMatrix('R', 'X', 0.00462963);
  SetMatrix('R', '*', -0.106481);
  SetMatrix('N', 'A', -0.0583333);
  SetMatrix('N', 'R', 0.0416667);
  SetMatrix('N', 'N', 0.341667);
  SetMatrix('N', 'D', 0.0916667);
  SetMatrix('N', 'C', -0.108333);
  SetMatrix('N', 'Q', 0.0416667);
  SetMatrix('N', 'E', 0.0416667);
  SetMatrix('N', 'G', 0.0416667);
  SetMatrix('N', 'H', 0.0916667);
  SetMatrix('N', 'I', -0.108333);
  SetMatrix('N', 'L', -0.108333);
  SetMatrix('N', 'K', 0.0416667);
  SetMatrix('N', 'M', -0.0583333);
  SetMatrix('N', 'F', -0.108333);
  SetMatrix('N', 'P', -0.0583333);
  SetMatrix('N', 'S', 0.0916667);
  SetMatrix('N', 'T', 0.0416667);
  SetMatrix('N', 'W', -0.158333);
  SetMatrix('N', 'Y', -0.0583333);
  SetMatrix('N', 'V', -0.108333);
  SetMatrix('N', 'B', 0.191667);
  SetMatrix('N', 'Z', 0.0416667);
  SetMatrix('N', 'X', -0.00833333);
  SetMatrix('N', '*', -0.158333);
  SetMatrix('D', 'A', -0.0352564);
  SetMatrix('D', 'R', -0.0352564);
  SetMatrix('D', 'N', 0.0801282);
  SetMatrix('D', 'D', 0.272436);
  SetMatrix('D', 'C', -0.0737179);
  SetMatrix('D', 'Q', 0.0416667);
  SetMatrix('D', 'E', 0.11859);
  SetMatrix('D', 'G', 0.00320513);
  SetMatrix('D', 'H', 0.00320513);
  SetMatrix('D', 'I', -0.0737179);
  SetMatrix('D', 'L', -0.112179);
  SetMatrix('D', 'K', 0.00320513);
  SetMatrix('D', 'M', -0.0737179);
  SetMatrix('D', 'F', -0.0737179);
  SetMatrix('D', 'P', 0.00320513);
  SetMatrix('D', 'S', 0.0416667);
  SetMatrix('D', 'T', 0.00320513);
  SetMatrix('D', 'W', -0.112179);
  SetMatrix('D', 'Y', -0.0737179);
  SetMatrix('D', 'V', -0.0737179);
  SetMatrix('D', 'B', 0.195513);
  SetMatrix('D', 'Z', 0.0801282);
  SetMatrix('D', 'X', 0.00320513);
  SetMatrix('D', '*', -0.112179);
  SetMatrix('C', 'A', 0.0416667);
  SetMatrix('C', 'R', -0.0281008);
  SetMatrix('C', 'N', -0.0281008);
  SetMatrix('C', 'D', -0.0281008);
  SetMatrix('C', 'C', 0.250969);
  SetMatrix('C', 'Q', -0.0281008);
  SetMatrix('C', 'E', -0.0513566);
  SetMatrix('C', 'G', -0.0281008);
  SetMatrix('C', 'H', -0.0281008);
  SetMatrix('C', 'I', 0.0184109);
  SetMatrix('C', 'L', 0.0184109);
  SetMatrix('C', 'K', -0.0281008);
  SetMatrix('C', 'M', 0.0184109);
  SetMatrix('C', 'F', -0.00484496);
  SetMatrix('C', 'P', -0.0281008);
  SetMatrix('C', 'S', 0.0184109);
  SetMatrix('C', 'T', 0.0184109);
  SetMatrix('C', 'W', -0.00484496);
  SetMatrix('C', 'Y', -0.00484496);
  SetMatrix('C', 'V', 0.0184109);
  SetMatrix('C', 'B', -0.0281008);
  SetMatrix('C', 'Z', -0.0281008);
  SetMatrix('C', 'X', -0.00484496);
  SetMatrix('C', '*', -0.0513566);
  SetMatrix('Q', 'A', -0.0297619);
  SetMatrix('Q', 'R', 0.113095);
  SetMatrix('Q', 'N', 0.0416667);
  SetMatrix('Q', 'D', 0.0416667);
  SetMatrix('Q', 'C', -0.172619);
  SetMatrix('Q', 'Q', 0.39881);
  SetMatrix('Q', 'E', 0.184524);
  SetMatrix('Q', 'G', -0.10119);
  SetMatrix('Q', 'H', 0.0416667);
  SetMatrix('Q', 'I', -0.172619);
  SetMatrix('Q', 'L', -0.10119);
  SetMatrix('Q', 'K', 0.113095);
  SetMatrix('Q', 'M', 0.0416667);
  SetMatrix('Q', 'F', -0.172619);
  SetMatrix('Q', 'P', -0.0297619);
  SetMatrix('Q', 'S', 0.0416667);
  SetMatrix('Q', 'T', -0.0297619);
  SetMatrix('Q', 'W', -0.10119);
  SetMatrix('Q', 'Y', -0.0297619);
  SetMatrix('Q', 'V', -0.10119);
  SetMatrix('Q', 'B', 0.0416667);
  SetMatrix('Q', 'Z', 0.255952);
  SetMatrix('Q', 'X', -0.0297619);
  SetMatrix('Q', '*', -0.244048);
  SetMatrix('E', 'A', -0.0171569);
  SetMatrix('E', 'R', 0.0416667);
  SetMatrix('E', 'N', 0.0416667);
  SetMatrix('E', 'D', 0.159314);
  SetMatrix('E', 'C', -0.193627);
  SetMatrix('E', 'Q', 0.159314);
  SetMatrix('E', 'E', 0.335784);
  SetMatrix('E', 'G', -0.0759804);
  SetMatrix('E', 'H', 0.0416667);
  SetMatrix('E', 'I', -0.134804);
  SetMatrix('E', 'L', -0.134804);
  SetMatrix('E', 'K', 0.10049);
  SetMatrix('E', 'M', -0.0759804);
  SetMatrix('E', 'F', -0.134804);
  SetMatrix('E', 'P', -0.0171569);
  SetMatrix('E', 'S', 0.0416667);
  SetMatrix('E', 'T', -0.0171569);
  SetMatrix('E', 'W', -0.134804);
  SetMatrix('E', 'Y', -0.0759804);
  SetMatrix('E', 'V', -0.0759804);
  SetMatrix('E', 'B', 0.10049);
  SetMatrix('E', 'Z', 0.276961);
  SetMatrix('E', 'X', -0.0171569);
  SetMatrix('E', '*', -0.193627);
  SetMatrix('G', 'A', 0.0416667);
  SetMatrix('G', 'R', -0.00595238);
  SetMatrix('G', 'N', 0.0416667);
  SetMatrix('G', 'D', 0.0178571);
  SetMatrix('G', 'C', -0.0297619);
  SetMatrix('G', 'Q', -0.00595238);
  SetMatrix('G', 'E', -0.00595238);
  SetMatrix('G', 'G', 0.184524);
  SetMatrix('G', 'H', -0.00595238);
  SetMatrix('G', 'I', -0.0535714);
  SetMatrix('G', 'L', -0.0535714);
  SetMatrix('G', 'K', -0.00595238);
  SetMatrix('G', 'M', -0.0297619);
  SetMatrix('G', 'F', -0.0297619);
  SetMatrix('G', 'P', -0.00595238);
  SetMatrix('G', 'S', 0.0416667);
  SetMatrix('G', 'T', -0.00595238);
  SetMatrix('G', 'W', -0.00595238);
  SetMatrix('G', 'Y', -0.0297619);
  SetMatrix('G', 'V', -0.0297619);
  SetMatrix('G', 'B', 0.0178571);
  SetMatrix('G', 'Z', -0.00595238);
  SetMatrix('G', 'X', 0.0178571);
  SetMatrix('G', '*', -0.0535714);
  SetMatrix('H', 'A', -0.0492424);
  SetMatrix('H', 'R', 0.0416667);
  SetMatrix('H', 'N', 0.0871212);
  SetMatrix('H', 'D', -0.00378788);
  SetMatrix('H', 'C', -0.094697);
  SetMatrix('H', 'Q', 0.0416667);
  SetMatrix('H', 'E', 0.0416667);
  SetMatrix('H', 'G', -0.0492424);
  SetMatrix('H', 'H', 0.405303);
  SetMatrix('H', 'I', -0.094697);
  SetMatrix('H', 'L', -0.094697);
  SetMatrix('H', 'K', -0.00378788);
  SetMatrix('H', 'M', -0.0492424);
  SetMatrix('H', 'F', -0.00378788);
  SetMatrix('H', 'P', -0.0492424);
  SetMatrix('H', 'S', -0.00378788);
  SetMatrix('H', 'T', -0.0492424);
  SetMatrix('H', 'W', -0.0492424);
  SetMatrix('H', 'Y', 0.132576);
  SetMatrix('H', 'V', -0.094697);
  SetMatrix('H', 'B', 0.0416667);
  SetMatrix('H', 'Z', 0.0416667);
  SetMatrix('H', 'X', -0.00378788);
  SetMatrix('H', '*', -0.140152);
  SetMatrix('I', 'A', 0.0153509);
  SetMatrix('I', 'R', -0.0372807);
  SetMatrix('I', 'N', -0.0372807);
  SetMatrix('I', 'D', -0.0372807);
  SetMatrix('I', 'C', 0.0153509);
  SetMatrix('I', 'Q', -0.0372807);
  SetMatrix('I', 'E', -0.0372807);
  SetMatrix('I', 'G', -0.0635965);
  SetMatrix('I', 'H', -0.0372807);
  SetMatrix('I', 'I', 0.14693);
  SetMatrix('I', 'L', 0.0942982);
  SetMatrix('I', 'K', -0.0372807);
  SetMatrix('I', 'M', 0.0679825);
  SetMatrix('I', 'F', 0.0416667);
  SetMatrix('I', 'P', -0.0372807);
  SetMatrix('I', 'S', -0.0109649);
  SetMatrix('I', 'T', 0.0153509);
  SetMatrix('I', 'W', -0.0372807);
  SetMatrix('I', 'Y', 0.0153509);
  SetMatrix('I', 'V', 0.120614);
  SetMatrix('I', 'B', -0.0372807);
  SetMatrix('I', 'Z', -0.0372807);
  SetMatrix('I', 'X', 0.0153509);
  SetMatrix('I', '*', -0.0635965);
  SetMatrix('L', 'A', 0.0146396);
  SetMatrix('L', 'R', -0.0123874);
  SetMatrix('L', 'N', -0.0394144);
  SetMatrix('L', 'D', -0.0664414);
  SetMatrix('L', 'C', 0.0146396);
  SetMatrix('L', 'Q', -0.0123874);
  SetMatrix('L', 'E', -0.0394144);
  SetMatrix('L', 'G', -0.0664414);
  SetMatrix('L', 'H', -0.0394144);
  SetMatrix('L', 'I', 0.0957207);
  SetMatrix('L', 'L', 0.149775);
  SetMatrix('L', 'K', -0.0123874);
  SetMatrix('L', 'M', 0.0957207);
  SetMatrix('L', 'F', 0.0416667);
  SetMatrix('L', 'P', -0.0394144);
  SetMatrix('L', 'S', -0.0123874);
  SetMatrix('L', 'T', 0.0146396);
  SetMatrix('L', 'W', -0.0123874);
  SetMatrix('L', 'Y', 0.0146396);
  SetMatrix('L', 'V', 0.0686937);
  SetMatrix('L', 'B', -0.0664414);
  SetMatrix('L', 'Z', -0.0394144);
  SetMatrix('L', 'X', 0.0146396);
  SetMatrix('L', '*', -0.0664414);
  SetMatrix('K', 'A', -0.00595238);
  SetMatrix('K', 'R', 0.136905);
  SetMatrix('K', 'N', 0.0416667);
  SetMatrix('K', 'D', -0.00595238);
  SetMatrix('K', 'C', -0.10119);
  SetMatrix('K', 'Q', 0.0892857);
  SetMatrix('K', 'E', 0.0892857);
  SetMatrix('K', 'G', -0.0535714);
  SetMatrix('K', 'H', -0.00595238);
  SetMatrix('K', 'I', -0.10119);
  SetMatrix('K', 'L', -0.0535714);
  SetMatrix('K', 'K', 0.279762);
  SetMatrix('K', 'M', -0.00595238);
  SetMatrix('K', 'F', -0.10119);
  SetMatrix('K', 'P', -0.00595238);
  SetMatrix('K', 'S', 0.0416667);
  SetMatrix('K', 'T', -0.00595238);
  SetMatrix('K', 'W', -0.10119);
  SetMatrix('K', 'Y', -0.0535714);
  SetMatrix('K', 'V', -0.0535714);
  SetMatrix('K', 'B', 0.0416667);
  SetMatrix('K', 'Z', 0.0892857);
  SetMatrix('K', 'X', -0.00595238);
  SetMatrix('K', '*', -0.14881);
  SetMatrix('M', 'A', -0.00378788);
  SetMatrix('M', 'R', -0.00378788);
  SetMatrix('M', 'N', -0.0492424);
  SetMatrix('M', 'D', -0.094697);
  SetMatrix('M', 'C', -0.00378788);
  SetMatrix('M', 'Q', 0.0416667);
  SetMatrix('M', 'E', -0.0492424);
  SetMatrix('M', 'G', -0.094697);
  SetMatrix('M', 'H', -0.0492424);
  SetMatrix('M', 'I', 0.0871212);
  SetMatrix('M', 'L', 0.132576);
  SetMatrix('M', 'K', -0.00378788);
  SetMatrix('M', 'M', 0.268939);
  SetMatrix('M', 'F', 0.0416667);
  SetMatrix('M', 'P', -0.0492424);
  SetMatrix('M', 'S', -0.00378788);
  SetMatrix('M', 'T', -0.00378788);
  SetMatrix('M', 'W', -0.00378788);
  SetMatrix('M', 'Y', -0.00378788);
  SetMatrix('M', 'V', 0.0871212);
  SetMatrix('M', 'B', -0.094697);
  SetMatrix('M', 'Z', -0.00378788);
  SetMatrix('M', 'X', -0.00378788);
  SetMatrix('M', '*', -0.140152);
  SetMatrix('F', 'A', -0.0138889);
  SetMatrix('F', 'R', -0.0416667);
  SetMatrix('F', 'N', -0.0416667);
  SetMatrix('F', 'D', -0.0416667);
  SetMatrix('F', 'C', -0.0138889);
  SetMatrix('F', 'Q', -0.0416667);
  SetMatrix('F', 'E', -0.0416667);
  SetMatrix('F', 'G', -0.0416667);
  SetMatrix('F', 'H', 0.0138889);
  SetMatrix('F', 'I', 0.0416667);
  SetMatrix('F', 'L', 0.0416667);
  SetMatrix('F', 'K', -0.0416667);
  SetMatrix('F', 'M', 0.0416667);
  SetMatrix('F', 'F', 0.208333);
  SetMatrix('F', 'P', -0.0694444);
  SetMatrix('F', 'S', -0.0138889);
  SetMatrix('F', 'T', -0.0138889);
  SetMatrix('F', 'W', 0.0694444);
  SetMatrix('F', 'Y', 0.125);
  SetMatrix('F', 'V', 0.0138889);
  SetMatrix('F', 'B', -0.0416667);
  SetMatrix('F', 'Z', -0.0416667);
  SetMatrix('F', 'X', 0.0138889);
  SetMatrix('F', '*', -0.0694444);
  SetMatrix('P', 'A', 0.0172764);
  SetMatrix('P', 'R', -0.00711382);
  SetMatrix('P', 'N', -0.00711382);
  SetMatrix('P', 'D', 0.0172764);
  SetMatrix('P', 'C', -0.0315041);
  SetMatrix('P', 'Q', 0.0172764);
  SetMatrix('P', 'E', 0.0172764);
  SetMatrix('P', 'G', -0.00711382);
  SetMatrix('P', 'H', -0.00711382);
  SetMatrix('P', 'I', -0.0315041);
  SetMatrix('P', 'L', -0.0315041);
  SetMatrix('P', 'K', 0.0172764);
  SetMatrix('P', 'M', -0.00711382);
  SetMatrix('P', 'F', -0.0558943);
  SetMatrix('P', 'P', 0.212398);
  SetMatrix('P', 'S', 0.0172764);
  SetMatrix('P', 'T', 0.0172764);
  SetMatrix('P', 'W', -0.0558943);
  SetMatrix('P', 'Y', -0.0315041);
  SetMatrix('P', 'V', -0.00711382);
  SetMatrix('P', 'B', -0.00711382);
  SetMatrix('P', 'Z', 0.0172764);
  SetMatrix('P', 'X', -0.00711382);
  SetMatrix('P', '*', -0.0558943);
  SetMatrix('S', 'A', 0.108333);
  SetMatrix('S', 'R', -0.025);
  SetMatrix('S', 'N', 0.108333);
  SetMatrix('S', 'D', 0.0416667);
  SetMatrix('S', 'C', -0.025);
  SetMatrix('S', 'Q', 0.0416667);
  SetMatrix('S', 'E', 0.0416667);
  SetMatrix('S', 'G', 0.0416667);
  SetMatrix('S', 'H', -0.025);
  SetMatrix('S', 'I', -0.0916667);
  SetMatrix('S', 'L', -0.0916667);
  SetMatrix('S', 'K', 0.0416667);
  SetMatrix('S', 'M', -0.025);
  SetMatrix('S', 'F', -0.0916667);
  SetMatrix('S', 'P', -0.025);
  SetMatrix('S', 'S', 0.308333);
  SetMatrix('S', 'T', 0.108333);
  SetMatrix('S', 'W', -0.158333);
  SetMatrix('S', 'Y', -0.0916667);
  SetMatrix('S', 'V', -0.0916667);
  SetMatrix('S', 'B', 0.0416667);
  SetMatrix('S', 'Z', 0.0416667);
  SetMatrix('S', 'X', 0.0416667);
  SetMatrix('S', '*', -0.225);
  SetMatrix('T', 'A', 0.0416667);
  SetMatrix('T', 'R', -0.00833333);
  SetMatrix('T', 'N', 0.0416667);
  SetMatrix('T', 'D', -0.00833333);
  SetMatrix('T', 'C', -0.00833333);
  SetMatrix('T', 'Q', -0.00833333);
  SetMatrix('T', 'E', -0.00833333);
  SetMatrix('T', 'G', -0.0583333);
  SetMatrix('T', 'H', -0.0583333);
  SetMatrix('T', 'I', -0.00833333);
  SetMatrix('T', 'L', -0.00833333);
  SetMatrix('T', 'K', -0.00833333);
  SetMatrix('T', 'M', -0.00833333);
  SetMatrix('T', 'F', -0.0583333);
  SetMatrix('T', 'P', -0.00833333);
  SetMatrix('T', 'S', 0.0916667);
  SetMatrix('T', 'T', 0.291667);
  SetMatrix('T', 'W', -0.0583333);
  SetMatrix('T', 'Y', -0.0583333);
  SetMatrix('T', 'V', 0.0416667);
  SetMatrix('T', 'B', -0.00833333);
  SetMatrix('T', 'Z', -0.00833333);
  SetMatrix('T', 'X', 0.0416667);
  SetMatrix('T', '*', -0.158333);
  SetMatrix('W', 'A', -0.025);
  SetMatrix('W', 'R', -0.025);
  SetMatrix('W', 'N', -0.0472222);
  SetMatrix('W', 'D', -0.0472222);
  SetMatrix('W', 'C', -0.00277778);
  SetMatrix('W', 'Q', -0.00277778);
  SetMatrix('W', 'E', -0.025);
  SetMatrix('W', 'G', -0.00277778);
  SetMatrix('W', 'H', -0.00277778);
  SetMatrix('W', 'I', -0.025);
  SetMatrix('W', 'L', -0.00277778);
  SetMatrix('W', 'K', -0.025);
  SetMatrix('W', 'M', 0.0194444);
  SetMatrix('W', 'F', 0.0638889);
  SetMatrix('W', 'P', -0.0472222);
  SetMatrix('W', 'S', -0.025);
  SetMatrix('W', 'T', -0.00277778);
  SetMatrix('W', 'W', 0.286111);
  SetMatrix('W', 'Y', 0.0861111);
  SetMatrix('W', 'V', -0.025);
  SetMatrix('W', 'B', -0.0472222);
  SetMatrix('W', 'Z', -0.025);
  SetMatrix('W', 'X', -0.00277778);
  SetMatrix('W', '*', -0.0472222);
  SetMatrix('Y', 'A', -0.0352564);
  SetMatrix('Y', 'R', -0.0352564);
  SetMatrix('Y', 'N', -0.0352564);
  SetMatrix('Y', 'D', -0.0737179);
  SetMatrix('Y', 'C', -0.0352564);
  SetMatrix('Y', 'Q', 0.00320513);
  SetMatrix('Y', 'E', -0.0352564);
  SetMatrix('Y', 'G', -0.0737179);
  SetMatrix('Y', 'H', 0.11859);
  SetMatrix('Y', 'I', 0.00320513);
  SetMatrix('Y', 'L', 0.00320513);
  SetMatrix('Y', 'K', -0.0352564);
  SetMatrix('Y', 'M', 0.00320513);
  SetMatrix('Y', 'F', 0.157051);
  SetMatrix('Y', 'P', -0.0737179);
  SetMatrix('Y', 'S', -0.0352564);
  SetMatrix('Y', 'T', -0.0352564);
  SetMatrix('Y', 'W', 0.11859);
  SetMatrix('Y', 'Y', 0.310897);
  SetMatrix('Y', 'V', 0.00320513);
  SetMatrix('Y', 'B', -0.0737179);
  SetMatrix('Y', 'Z', -0.0352564);
  SetMatrix('Y', 'X', 0.00320513);
  SetMatrix('Y', '*', -0.112179);
  SetMatrix('V', 'A', 0.0416667);
  SetMatrix('V', 'R', -0.0520833);
  SetMatrix('V', 'N', -0.0520833);
  SetMatrix('V', 'D', -0.0520833);
  SetMatrix('V', 'C', 0.0104167);
  SetMatrix('V', 'Q', -0.0208333);
  SetMatrix('V', 'E', -0.0208333);
  SetMatrix('V', 'G', -0.0520833);
  SetMatrix('V', 'H', -0.0520833);
  SetMatrix('V', 'I', 0.135417);
  SetMatrix('V', 'L', 0.0729167);
  SetMatrix('V', 'K', -0.0208333);
  SetMatrix('V', 'M', 0.0729167);
  SetMatrix('V', 'F', 0.0104167);
  SetMatrix('V', 'P', -0.0208333);
  SetMatrix('V', 'S', -0.0208333);
  SetMatrix('V', 'T', 0.0416667);
  SetMatrix('V', 'W', -0.0520833);
  SetMatrix('V', 'Y', 0.0104167);
  SetMatrix('V', 'V', 0.166667);
  SetMatrix('V', 'B', -0.0520833);
  SetMatrix('V', 'Z', -0.0208333);
  SetMatrix('V', 'X', 0.0104167);
  SetMatrix('V', '*', -0.0833333);
  SetMatrix('B', 'A', -0.0383333);
  SetMatrix('B', 'R', 0.00166667);
  SetMatrix('B', 'N', 0.161667);
  SetMatrix('B', 'D', 0.201667);
  SetMatrix('B', 'C', -0.0783333);
  SetMatrix('B', 'Q', 0.0416667);
  SetMatrix('B', 'E', 0.0816667);
  SetMatrix('B', 'G', 0.00166667);
  SetMatrix('B', 'H', 0.0416667);
  SetMatrix('B', 'I', -0.0783333);
  SetMatrix('B', 'L', -0.118333);
  SetMatrix('B', 'K', 0.0416667);
  SetMatrix('B', 'M', -0.0783333);
  SetMatrix('B', 'F', -0.0783333);
  SetMatrix('B', 'P', -0.0383333);
  SetMatrix('B', 'S', 0.0416667);
  SetMatrix('B', 'T', 0.00166667);
  SetMatrix('B', 'W', -0.118333);
  SetMatrix('B', 'Y', -0.0783333);
  SetMatrix('B', 'V', -0.0783333);
  SetMatrix('B', 'B', 0.201667);
  SetMatrix('B', 'Z', 0.0816667);
  SetMatrix('B', 'X', 0.00166667);
  SetMatrix('B', '*', -0.118333);
  SetMatrix('Z', 'A', -0.0208333);
  SetMatrix('Z', 'R', 0.0416667);
  SetMatrix('Z', 'N', 0.0416667);
  SetMatrix('Z', 'D', 0.104167);
  SetMatrix('Z', 'C', -0.145833);
  SetMatrix('Z', 'Q', 0.229167);
  SetMatrix('Z', 'E', 0.291667);
  SetMatrix('Z', 'G', -0.0833333);
  SetMatrix('Z', 'H', 0.0416667);
  SetMatrix('Z', 'I', -0.145833);
  SetMatrix('Z', 'L', -0.145833);
  SetMatrix('Z', 'K', 0.104167);
  SetMatrix('Z', 'M', -0.0208333);
  SetMatrix('Z', 'F', -0.145833);
  SetMatrix('Z', 'P', -0.0208333);
  SetMatrix('Z', 'S', 0.0416667);
  SetMatrix('Z', 'T', -0.0208333);
  SetMatrix('Z', 'W', -0.145833);
  SetMatrix('Z', 'Y', -0.0833333);
  SetMatrix('Z', 'V', -0.0833333);
  SetMatrix('Z', 'B', 0.104167);
  SetMatrix('Z', 'Z', 0.291667);
  SetMatrix('Z', 'X', -0.0208333);
  SetMatrix('Z', '*', -0.208333);
  SetMatrix('X', 'A', 0.0416667);
  SetMatrix('X', 'R', 0.00462963);
  SetMatrix('X', 'N', 0.00462963);
  SetMatrix('X', 'D', 0.00462963);
  SetMatrix('X', 'C', -0.0324074);
  SetMatrix('X', 'Q', 0.00462963);
  SetMatrix('X', 'E', 0.00462963);
  SetMatrix('X', 'G', 0.00462963);
  SetMatrix('X', 'H', 0.00462963);
  SetMatrix('X', 'I', 0.00462963);
  SetMatrix('X', 'L', 0.00462963);
  SetMatrix('X', 'K', 0.00462963);
  SetMatrix('X', 'M', 0.00462963);
  SetMatrix('X', 'F', 0.00462963);
  SetMatrix('X', 'P', -0.0324074);
  SetMatrix('X', 'S', 0.0416667);
  SetMatrix('X', 'T', 0.0416667);
  SetMatrix('X', 'W', -0.0324074);
  SetMatrix('X', 'Y', 0.00462963);
  SetMatrix('X', 'V', 0.00462963);
  SetMatrix('X', 'B', 0.00462963);
  SetMatrix('X', 'Z', 0.00462963);
  SetMatrix('X', 'X', 0.00462963);
  SetMatrix('X', '*', -0.106481);
  SetMatrix('*', 'A', -0.00228938);
  SetMatrix('*', 'R', -0.00228938);
  SetMatrix('*', 'N', -0.00228938);
  SetMatrix('*', 'D', -0.00228938);
  SetMatrix('*', 'C', -0.00228938);
  SetMatrix('*', 'Q', -0.00228938);
  SetMatrix('*', 'E', -0.00228938);
  SetMatrix('*', 'G', -0.00228938);
  SetMatrix('*', 'H', -0.00228938);
  SetMatrix('*', 'I', -0.00228938);
  SetMatrix('*', 'L', -0.00228938);
  SetMatrix('*', 'K', -0.00228938);
  SetMatrix('*', 'M', -0.00228938);
  SetMatrix('*', 'F', -0.00228938);
  SetMatrix('*', 'P', -0.00228938);
  SetMatrix('*', 'S', -0.00228938);
  SetMatrix('*', 'T', -0.00228938);
  SetMatrix('*', 'W', -0.00228938);
  SetMatrix('*', 'Y', -0.00228938);
  SetMatrix('*', 'V', -0.00228938);
  SetMatrix('*', 'B', -0.00228938);
  SetMatrix('*', 'Z', -0.00228938);
  SetMatrix('*', 'X', -0.00228938);
  SetMatrix('*', '*', 0.0526557);



}
