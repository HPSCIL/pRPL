#ifndef  PGTIOL_GTIFFMETADATA_H
#define  PGTIOL_GTIFFMETADATA_H 

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <fcntl.h> 
#include <sstream>
#include <iostream>
#include <string>
#include <gdal_priv.h>
#include <cpl_string.h>
#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <tiff.h>
#include <tiffio.h>
#include <mpi.h>
#include "pgtiol-basicTypes.h"

using namespace std;

class gTiffMetaData {
  public:
    gTiffMetaData(void);
    ~gTiffMetaData(void);

   
    geoTiffMetaData* getMetaData() {return tiffMetaData;}
    void setProjection(string projection) {tiffMetaData->projectionRef = projection;}
    const char* getProjectionRef(void) {return tiffMetaData->projectionRef.c_str();}
    uint32_t getRasterYSize() {return tiffMetaData->y_size;}
    uint32_t getRasterXSize() {return tiffMetaData->x_size;}
    double getNoDataValue() {return tiffMetaData->no_data;}
    void setNoDataValue(double nodata) {tiffMetaData->no_data = nodata;}

    long getTileWidth() {return (long)tiffMetaData->blocks_x_size;}
    long getTileLength() {return (long)tiffMetaData->blocks_y_size;}
    void setTileWidth(long tileWidth) { tiffMetaData->blocks_x_size = tileWidth;}
    void setTileLength(long tileLength) { tiffMetaData->blocks_y_size = tileLength;}
    unsigned short getSampleFormat() {return tiffMetaData->sampleFormat;}
    unsigned short getDataSizeFileIn() {return tiffMetaData->dataSizeFileIn;}

    void getGeoTransform(double aGeoTransform[6]) {
      aGeoTransform[0] = tiffMetaData->geoTransform[0];
      aGeoTransform[1] = tiffMetaData->geoTransform[1];
      aGeoTransform[2] = tiffMetaData->geoTransform[2];
      aGeoTransform[3] = tiffMetaData->geoTransform[3];
      aGeoTransform[4] = tiffMetaData->geoTransform[4];
      aGeoTransform[5] = tiffMetaData->geoTransform[5];
    }
    void setGeoTransform(double aGeoTransform[6]){
      tiffMetaData->geoTransform[0] =aGeoTransform[0];
      tiffMetaData->geoTransform[1] =aGeoTransform[1];
      tiffMetaData->geoTransform[2] =aGeoTransform[2];
      tiffMetaData->geoTransform[3] =aGeoTransform[3];
      tiffMetaData->geoTransform[4] =aGeoTransform[4];
      tiffMetaData->geoTransform[5] =aGeoTransform[5];
    }
    void setTotalX(uint32_t totalx) {tiffMetaData->x_size = totalx;}
    void setTotalY(uint32_t totaly) {tiffMetaData->y_size = totaly;}
    int getBands() {return tiffMetaData->band_count;}
    void setBands(const int bands) {tiffMetaData->band_count = bands;}
    void setDataType(GDALDataType datatype) {
      tiffMetaData->band_type = datatype;
      tiffMetaData->dataSizeObj = getDataObj(datatype);
    }
    GDALDataType  getDataType(){return tiffMetaData->band_type;}

    bool setInputFile(char*fileName);
    bool readMetaData();

    bool setOutputFile(char*outputFile);
    bool writeMetaData();
   
    bool getMetaDataRaster();

  private:
    geoTiffMetaData *tiffMetaData;    
    char inputGTiffFile[MAXLEN];    
    char outputGTiffFile[MAXLEN];  

   
    int16_t parse_int16(uint8_t *buffer, bool big_endian);
    int64_t read_int64(MPI_File&fh, int64_t offset, bool big_endian);
    int64_t parse_int64(uint8_t *buffer, bool big_endian);
    string intToString(int output_tile_size);
    int write_int64(MPI_File&fh,int64_t offset,int64_t value,bool big_endian);
    int export_int64(int64_t num, uint8_t *buffer, bool big_endian);
    int getDataObj(GDALDataType band_type);
};

#endif

