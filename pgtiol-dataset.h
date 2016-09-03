#ifndef PGTIOL_DATASET_H
#define PGTIOL_DATASET_H
#include <gdal_priv.h>
#include <gdal.h>
#include "pgtiol-basicTypes.h"
#include "pgtiol-gTiffManager.h"

typedef enum {
  PG_ReadOnly = 0,
  PG_Update = 1
} PGTIOLAccess;

typedef enum {
  PG_None = 0,
  PG_Debug = 1,
  PG_Warning = 2,
  PG_Failure = 3,
  PG_Fatal = 4
} PGTIOLCPLErr;

typedef enum {
  PG_Read = 0,
  PG_Write = 1
} PGTIOLRWFlag;

class PGTIOLDriver;

class PGTIOLDataset {
public:
  PGTIOLDataset(const char *aFileName,
                PGTIOLAccess ioOption,
                char **papszOptions,
                int nXSize,
                int nYSize,
                int nBands,
                GDALDataType eType);
  PGTIOLDataset(const char *aFileName,
                PGTIOLAccess ioOption = PG_Update,
                int proc=0);
  ~PGTIOLDataset();

  int    GetRasterBand() { return m_tiffManager->getRasterCount();}
  void   SetRasterBand(const int bands) {m_tiffManager->setRasterCount(bands);}
  long   GetTileWidth(){return m_tiffManager->getTileWidth();}
  void   SetTileWidth(const long tileWidth){m_tiffManager->setTileWidth(tileWidth);}
  long   GetTileLength(){return m_tiffManager->getTileLength();}
  void   SetTileLength(const long tileLength){m_tiffManager->setTileLength(tileLength);}
  
  
  void SetProjection(string projection){ m_tiffManager->setProjection(projection);}
  void SetGeoTransform(double *GeoTransform){  m_tiffManager->setGeoTransform(GeoTransform);}
  void GetGeoTransform(double *GeoTransform){m_tiffManager->getGeoTransform(GeoTransform);}
  const char *GetProjectionRef(){return m_tiffManager->getProjectionRef();}
  void SetOutputGTiff(){ m_tiffManager->writeMetaData();}
 
  GDALDataType GetRasterDataType(){
  return m_tiffManager->getDataType();}
  int GetRasterXSize(){ return m_tiffManager->getRasterXSize();}
  int GetRasterYSize(){return m_tiffManager->getRasterYSize();}
  int GetRasterCount(){return m_tiffManager->getRasterCount();}
  double GetNoDataValue(){return m_tiffManager->getNoDataValue();}
  void SetNoDataValue(double NoDataValue){m_tiffManager->setNoDataValue(NoDataValue);}
  void SetWriteProcess(int proc_id){m_tiffManager->setProcess(proc_id);}

  
  PGTIOLCPLErr RasterIO(PGTIOLRWFlag eRWFlag,int nXOff,int nYOff,int nXSize,int nYSize,void * pData,int nBufXSize,int nBufYSize,
  GDALDataType eBufType,int band_count,int nPixelSpace,int nLineSpace );
  bool Open(const char *aFileName);
private:
  GTiffManager *m_tiffManager;
};

class PGTIOLDriver {
public:
  PGTIOLDriver(){};
  ~PGTIOLDriver(){};
  PGTIOLDataset *Create(const char * pszName,
                        int nXSize, int nYSize, int nBands,
                        GDALDataType eType,
                        char **papszOptions );
};

static void PGTIOLClose(PGTIOLDataset *pt){
	delete pt;
}

static PGTIOLDataset* PGTIOLOpen(char *aFileName,
                                 PGTIOLAccess ioOption){
   
   return new PGTIOLDataset(aFileName,ioOption);
}

#endif
