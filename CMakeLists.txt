cmake_minimum_required(VERSION 3.3)
project(cola)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14 -O3 -w")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)

# include directory in find path where all dependency modules exist
include_directories(./)


# cola binaries
set(SOURCE_FILES_MULTIALIGNCOLA  ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc src/cola/AlignmentCola.cc src/cola/AlignmentSet.cc src/cola/Cola.cc src/cola/EditGraph.cc src/cola/MultiAlignCola.cc src/cola/NOIAligner.cc src/cola/NSGAaligner.cc src/cola/NSaligner.cc src/cola/SWGAaligner.cc ryggrad/src/util/mutil.cc) 
#set(SOURCE_FILES_FINDCOLASIGDIST  ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc src/cola/AlignmentCola.cc src/cola/findColaSigDist.cc ryggrad/src/util/mutil.cc) 
set(SOURCE_FILES_RUNCOLA  ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc src/cola/AlignmentCola.cc src/cola/AlignmentSet.cc src/cola/Cola.cc src/cola/EditGraph.cc src/cola/NSGAaligner.cc src/cola/NSaligner.cc src/cola/SWGAaligner.cc src/cola/runCola.cc ryggrad/src/util/mutil.cc) 
set(SOURCE_FILES_RUNFALIGN  ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc src/cola/AlignmentCola.cc src/cola/AlignmentSet.cc src/cola/Cola.cc src/cola/EditGraph.cc src/cola/NOIAligner.cc src/cola/NSGAaligner.cc src/cola/NSaligner.cc src/cola/SWGAaligner.cc src/fastAlign/AlignmentThreads.cc src/fastAlign/DNASeqs.cc src/fastAlign/FastAlignUnit.cc src/fastAlign/SeedingThreads.cc src/fastAlign/runFAlign.cc src/fastAlign/SeedingObjects.cc src/fastAlign/SyntenicSeeds.cc ryggrad/src/util/mutil.cc) 

add_executable(MultiAlignCola  ${SOURCE_FILES_MULTIALIGNCOLA})
#add_executable(FindColaSigDist ${SOURCE_FILES_FINDCOLASIGDIST})
add_executable(RunCola ${SOURCE_FILES_RUNCOLA})
add_executable(RunFAlign ${SOURCE_FILES_RUNFALIGN})






