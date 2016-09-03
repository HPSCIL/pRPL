#include "pgtiol-dataset.h"
#include <assert.h>

PGTIOLDataset::
PGTIOLDataset(const char *aFileName,
              PGTIOLAccess ioOption,
              char **papszOptions,
              int nXSize, int nYSize, int nBands,
              GDALDataType eType) {
  if(ioOption == PG_Update){
     m_tiffManager=new GTiffManager();
     m_tiffManager->setDataType(eType);
     m_tiffManager->setTotalX(nXSize);
     m_tiffManager->setTotalY(nYSize);
	 m_tiffManager->setRasterCount(nBands);
	 m_tiffManager->setOutputFile(const_cast<char*>(aFileName));
  }
  else{
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to open the geotiff file";
  }
}

PGTIOLDataset::
PGTIOLDataset(const char *aFileName,
              PGTIOLAccess ioOption,
              int proc){
  if(ioOption == PG_ReadOnly){
    m_tiffManager=new GTiffManager();
    //m_tiffManager->setProcess(proc);
    m_tiffManager->setInputFile(aFileName);
    if(!m_tiffManager->readMetaData()) {
      delete m_tiffManager;
      m_tiffManager = NULL;
    }
  }

 if(ioOption == PG_Update){
   m_tiffManager=new GTiffManager();
   m_tiffManager->setProcess(proc);
   m_tiffManager->setOutputFile(const_cast<char*>(aFileName));
 }
}


PGTIOLDataset::
~PGTIOLDataset() {
  if(m_tiffManager!=NULL)
     delete m_tiffManager;
}



PGTIOLCPLErr PGTIOLDataset::RasterIO(PGTIOLRWFlag eRWFlag,int nXOff,int nYOff,int nXSize,int nYSize,void * pData,int nBufXSize,int nBufYSize,
                                     GDALDataType eBufType,int band_count,int nPixelSpace,int nLineSpace){   
  if(eRWFlag==PG_Read){
	  if(nXOff>=0 && (nXOff+nXSize<=m_tiffManager->getRasterXSize()) && nXOff>=0 && (nYOff+nYSize<=m_tiffManager->getRasterYSize())){
	      m_tiffManager->readIO(nXOff, nYOff, nYSize,nXSize, pData,eBufType,band_count);
	}
    else
      return PG_Failure;
  }
  if(eRWFlag==PG_Write) {
      m_tiffManager->writeIO(nXOff, nYOff, nYSize, nXSize, pData,eBufType,nLineSpace);
  } 
  return PG_None;
}

   
PGTIOLDataset * PGTIOLDriver::Create( const char * pszName,
    int nXSize, int nYSize, int nBands,
    GDALDataType eType, char ** papszOptions ){
  
  return new PGTIOLDataset(const_cast<char*>(pszName),PG_Update, papszOptions,nXSize,nYSize, nBands,eType);

}
