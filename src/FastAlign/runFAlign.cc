#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string>
#include <omp.h>
#include "base/CommandLineParser.h"
#include "extern/logger/log.h"
#include "src/FastAlign/FastAlignUnit.h"


int main(int argc,char** argv)
{

    commandArg<string> a1Cmmd("-i","FASTA file containing target sequences that are to be aligned agaiinst the query data");
    commandArg<string> a2Cmmd("-q","FASTA file containing query sequences");
    commandArg<string> bCmmd("-o","File to Output alignments", "alignments.out");
    commandArg<int>    cCmmd("-b","Subread block step", 10);
    commandArg<int>    dCmmd("-S","Seed size for choosing candidates", 15);
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

    string targetSeqFile    = P.GetStringValueFor(a1Cmmd);
    string querySeqFile    = P.GetStringValueFor(a2Cmmd);
    string outFile         = P.GetStringValueFor(bCmmd);
    int    readBlockSize   = P.GetIntValueFor(cCmmd);
    int    seedSize        = P.GetIntValueFor(dCmmd);
    double minIdent        = P.GetDoubleValueFor(eCmmd);
    int    alignBand       = P.GetIntValueFor(fCmmd);
    string applicationFile = P.GetStringValueFor(gCmmd);
    int    numOfThreads    = P.GetIntValueFor(threadCmmd);
    
    FILE* pFile               = fopen(applicationFile.c_str(), "w");
    Output2FILE::Stream()     = pFile;
    FILELog::ReportingLevel() = logINFO; 
#if defined(FORCE_DEBUG)
//    FILELog::ReportingLevel() = logDEBUG4; 
#endif
    
#if defined(OPEN_MP)
    omp_set_num_threads(numOfThreads); //The sort functions still use OMP
#endif
    
    AlignmentParams params(readBlockSize, seedSize,
                           minIdent, alignBand, 0.1);
    FastAlignUnit FAUnit(targetSeqFile, querySeqFile, params, numOfThreads);
    ofstream fOut;
    fOut.open(outFile.c_str());
    FAUnit.alignAllSeqs(fOut);
    return 0;
}



