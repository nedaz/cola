#ifndef _ALIGNMENT_PARAMS_H_
#define _ALIGNMENT_PARAMS_H_
//======================================================
class AlignmentParams 
{
public:
    AlignmentParams(int stepSize=10, int seedSize=15, 
                    float minIdent=0.7, int alignmentBound=3, float minSeedCover=0.25)
                   :m_suffixStep(stepSize), m_seedSize(seedSize),
                    m_minIdent(minIdent), m_alignmentBound(alignmentBound), 
                    m_minSeedCover(minSeedCover)  { }

    int   getSuffixStep() const     { return m_suffixStep;     }  
    int   getSeedSize()  const      { return m_seedSize;       } 
    float getMinIdentity() const    { return m_minIdent;       }
    int   getAlignBand() const      { return m_alignmentBound; }
    float getMinSeedCover() const   { return m_minSeedCover;   }

    void  setSuffixStep(int sst)    { m_suffixStep     = sst;  }  
    void  setSeedSize(int ss)       { m_seedSize       = ss;   } 
    void  setMinIdentity(float idt) { m_minIdent       = idt;  }
    void  setAlignBand(int ab)      { m_alignmentBound = ab;   }
    void  setMinSeedCover(float sc) { m_minSeedCover   = sc;   }


private: 
    int     m_suffixStep;     /// Block step size used for constructing subreads
    int     m_seedSize;       /// Seed size for selecting candidate reads for assembly extension
    float   m_minIdent;       /// Minimum identity for accepting a candidate read as an assembly extension
    int     m_alignmentBound; /// Alignment bandwidth used for local alignment to decide on choosing candidate reads
    float   m_minSeedCover;   /// The minimum acceptance level of seed coverage for choosing alignment candidates 
};
//======================================================

#endif // _ALIGNMENT_PARAMS_H_
