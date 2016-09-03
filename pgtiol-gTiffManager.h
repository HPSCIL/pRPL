#ifndef  PGTIOL_GTIFFMANAGER_H
#define  PGTIOL_GTIFFMANAGER_H

#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <mpi.h>

#include  "pgtiol-basicTypes.h"
#include  "pgtiol-gTiffMetaData.h"
#include  "pgtiol-gTiffData.h"

using namespace std;

class GTiffManager {
  public:
    GTiffManager(void);
    ~GTiffManager(void);

    bool setInputFile(const char *inputFileName);
    bool setOutputFile(const char *outputFileName);

    uint32_t getRasterCount() {return igTiffMetaData->getBands();}
    void setRasterCount(const int bands) {return igTiffMetaData->setBands(bands);}

    void setProjection(string projection) {igTiffMetaData->setProjection(projection);}
    const char* getProjectionRef(void) {return igTiffMetaData->getProjectionRef();}

    uint32_t getRasterYSize() {return igTiffMetaData->getRasterYSize();}
    uint32_t getRasterXSize() {return igTiffMetaData->getRasterXSize();}
    void setTotalX(uint32_t totalx) {igTiffMetaData->setTotalX(totalx);}
    void setTotalY(uint32_t totaly) {igTiffMetaData->setTotalY(totaly);}

    double getNoDataValue() {return igTiffMetaData->getNoDataValue();}
    void setNoDataValue(double nodata) {igTiffMetaData->setNoDataValue(nodata);}

    void getGeoTransform(double aGeoTransform[6]) {igTiffMetaData->getGeoTransform(aGeoTransform);}
    void setGeoTransform(double aGeoTransform[6]) {igTiffMetaData->setGeoTransform(aGeoTransform);}

    long getTileWidth() {return igTiffMetaData->getTileWidth();}
    void setTileWidth(const long tileWidth) {igTiffMetaData->setTileWidth(tileWidth);}
    long getTileLength() {return igTiffMetaData->getTileLength();}
    void setTileLength(const long tileLength) {igTiffMetaData->setTileLength(tileLength);}

    GDALDataType getDataType() { return igTiffMetaData->getDataType();};
    void setDataType(GDALDataType datatype) {igTiffMetaData->setDataType(datatype);}

    bool readMetaData();
    void setProcess(const int procs) {writeProcs = procs;}
    bool writeMetaData();

    bool writeIO(long nXOff, long nYOff, long nYSize,long nXSize, void*pData,GDALDataType dataType,long nLineSpace);
    bool readIO(long nXOff, long nYOff, long nYSize,long nXSize, void*pData,GDALDataType dataType,int band_count);

  private:
    char inputFileName[MAXLEN]; 
    char outputFileName[MAXLEN];

    gTiffMetaData *igTiffMetaData;
    gTiffData *igTiffData;

    int size;
    int rank;
    int writeProcs;
};

#endif
