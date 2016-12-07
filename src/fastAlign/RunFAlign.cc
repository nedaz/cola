#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string>
#include <omp.h>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/Logger.h"
#include "FastAlignUnit.h"


int main(int argc,char** argv)
{

    commandArg<string> a1Cmmd("-q","FASTA file containing query sequences that are to be aligned against the reference target data");
    commandArg<string> a2Cmmd("-t","FASTA file containing the reference target sequences");
    commandArg<string> bCmmd("-o","File to Output alignments", "alignments.out");
    commandArg<int>    cCmmd("-b","Subread block step", 2);
    commandArg<int>    dCmmd("-S","Seed size for choosing candidates", 20);
    commandArg<double> eCmmd("-I","Minimum acceptable identity for seeding sequences", 0.4);
    commandArg<int>    fCmmd("-B","Bandwidth for local alignments", 3);
    commandArg<string> gCmmd("-L","Application logging file","application.log");
    commandArg<int>    threadCmmd("-T","Number of Cores to run with", 1);

    commandLineParser P(argc,argv);
    P.SetDescription("Seeding based aligner that filters sequences before running exhaustive alignment");
    
    P.registerArg(a1Cmmd);
    P.registerArg(a2Cmmd);
    P.registerArg(bCmmd);
    P.registerArg(cCmmd);
    P.registerArg(dCmmd);
    P.registerArg(eCmmd);
    P.registerArg(fCmmd);
    P.registerArg(gCmmd);
    P.registerArg(threadCmmd);
    P.parse();

    string querySeqFile   = P.GetStringValueFor(a1Cmmd);
    string targetSeqFile    = P.GetStringValueFor(a2Cmmd);
    string outFile         = P.GetStringValueFor(bCmmd);
    int    readBlockSize   = P.GetIntValueFor(cCmmd);
    int    seedSize        = P.GetIntValueFor(dCmmd);
    double minIdent        = P.GetDoubleValueFor(eCmmd);
    int    alignBand       = P.GetIntValueFor(fCmmd);
    string applicationFile = P.GetStringValueFor(gCmmd);
    int    numThreads      = P.GetIntValueFor(threadCmmd);
    
    FILE* pFile               = fopen(applicationFile.c_str(), "w");
    Output2FILE::Stream()     = pFile;
    FILELog::ReportingLevel() = logINFO; 
    
#if defined(OPEN_MP)
    omp_set_num_threads(numThreads); //The sort functions still use OMP
#endif
    
    ofstream fOut;
    fOut.open(outFile.c_str());

    AlignmentParams params(readBlockSize, seedSize,
                           minIdent, alignBand, 0.05); // TODO The seed coverage threshold needs to be looked into
    FastAlignTargetUnit qUnit(targetSeqFile, readBlockSize);

    FastAlignUnit FAUnit(querySeqFile, qUnit, params, numThreads);
    FAUnit.alignAllSeqs(numThreads, fOut);

    cout << "Reverse Complement of sequences:" << endl;
    FastAlignUnit revFAUnit(querySeqFile, qUnit, params, numThreads, true);
    revFAUnit.alignAllSeqs(numThreads, fOut);

    fOut.close();
    return 0;
}



