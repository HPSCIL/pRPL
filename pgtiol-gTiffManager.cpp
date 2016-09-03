#include "pgtiol-gTiffManager.h"

GTiffManager::GTiffManager(void)
{
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  igTiffMetaData = new gTiffMetaData();
  igTiffData = new gTiffData();
  writeProcs = -1;
}


GTiffManager::~GTiffManager(void)
{
  if(igTiffMetaData!=NULL){
    delete igTiffMetaData;
    igTiffMetaData = NULL;
  }
  if(igTiffData!=NULL){
    delete igTiffData;
    igTiffData = NULL;
  }
}


bool GTiffManager::setInputFile(const char *inputFile)
{
  if(inputFile==NULL){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: input filename is null (" <<inputFile<<")"<<endl;
    return false;
  }
  strcpy(inputFileName, inputFile);
  return true;
}

bool GTiffManager::setOutputFile(const char *outputFile)
{
  if(outputFile==NULL){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: output filename is null (" <<outputFile<<")"<<endl;
    return false;
  }
  strncpy(outputFileName,outputFile, MAXLEN);
  return true;
}

bool GTiffManager::readMetaData()
{   

  if(!igTiffMetaData->setInputFile(inputFileName)){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: input filename is null (" <<inputFileName<<")"<<endl;
    return false;
  }

  if(!igTiffMetaData->readMetaData()){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error occurred when reading file  (" <<inputFileName<<")"<<endl;
    return false;
  }

  return true;
}

bool GTiffManager::writeMetaData()
{
  if(!igTiffMetaData->setOutputFile(outputFileName)){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error output filename is null (" <<outputFileName<<")"<<endl;
    return false;
  }
  if(rank == writeProcs){
    if(!igTiffMetaData->writeMetaData()){
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error occurred when writing metadata of  (" <<outputFileName<<")"<<endl;
      return false;
    }
  }

  //MPI_Barrier(MPI_COMM_WORLD); // guan edit

  if(!igTiffMetaData->getMetaDataRaster()){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error occurred when reading metadata of (" <<outputFileName<<")"<<endl;
    return false;
  }
  return true;
}

bool GTiffManager::writeIO(long nXOff, long nYOff, long nYSize,long nXSize, void*pData,GDALDataType dataType,long nLineSpace)
{
  if(!igTiffData->setOutputFileName(outputFileName)){
    return false;
  }
  if(!igTiffData->writeData(nXOff, nYOff, nYSize, nXSize,pData,dataType,igTiffMetaData,nLineSpace)){
    return false;
  }
  return true;
}

bool GTiffManager::readIO(long nXOff, long nYOff, long nYSize,long nXSize, void*pData,GDALDataType dataType,int band_count)
{
  if(!igTiffData->setInputFileName(inputFileName)){
    return false;
  }
  if(!igTiffData->readData(nXOff, nYOff, nYSize, nXSize,pData,dataType,band_count,igTiffMetaData)){
    return false;
  }
  return true;
}

