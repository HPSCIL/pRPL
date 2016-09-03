#ifndef PRPL_LAYER_H
#define PRPL_LAYER_H

#include "prpl-basicTypes.h"
#include "prpl-cellspace.h"
#include "prpl-subCellspace.h"
#include "prpl-smplDcmp.h"

namespace pRPL {
  class Layer {
    public:
      /* Constructors and destructor */
      Layer(const char *aName = NULL);
      ~Layer();
      
      /* Information */
      const char* name() const;
      void name(const char *aName);
      const pRPL::SpaceDims* glbDims() const;
      const char* dataType() const;
      size_t dataSize() const;
      long tileSize() const;
	  
      const pRPL::CellspaceGeoinfo * glbGeoinfo() const;

      /* GDAL */
      bool openGdalDS(const char *aFileName,
                      GDALAccess ioOption = GA_ReadOnly);
      int dsBand() const;
      void dsBand(int iBand);
      void closeGdalDS();
      bool createGdalDS(const char *aFileName,
                        const char *aGdalFormat,
                        char **aGdalOptions = NULL);
      /*PGTIOL*/
      bool openPgtiolDS(const char *aFileName,
	                      int prc,
                        PGTIOLAccess ioOption = PG_ReadOnly);
      void closePgtiolDS();
      bool createPgtiolDS(const char *aFileName,
	                        int _prc,
                          char **aPgtOptions = NULL);
      
      /* Decompose & Distribution */
      bool decompose(int nRowSubspcs,
                     int nColSubspcs,
                     const pRPL::Neighborhood *pNbrhd = NULL);
      bool copyDcmp(const Layer *pFromLyr);
      bool isDecomposed() const;

      /* CellspaceInfo */
      void initCellspaceInfo();
      void initCellspaceInfo(const pRPL::SpaceDims &dims,
                             const char *aTypeName,
                             size_t typeSize,
                             const pRPL::CellspaceGeoinfo *pGeoinfo = NULL,
                             long  tileSize = TILEWIDTH);
      void initCellspaceInfo(const pRPL::CellspaceInfo &cellspcInfo);
      
	    bool initCellspaceInfoByGDAL();
	    bool initCellspaceInfoByPGTIOL();

      void delCellspaceInfo();
      pRPL::CellspaceInfo* cellspaceInfo();
      const pRPL::CellspaceInfo* cellspaceInfo() const;

      /* Cellspace */
      bool createCellspace(double *pInitVal = NULL);
	    bool loadCellspaceByGDAL();
      bool writeCellspaceByGDAL();
	    bool loadCellspaceByPGTIOL();
      bool writeCellspaceByPGTIOL();

      void delCellspace();
      pRPL::Cellspace* cellspace();
      const pRPL::Cellspace* cellspace() const;

      /* SubCellspaceInfo */
      int nSubCellspaceInfos() const;
      const list<pRPL::SubCellspaceInfo>* allSubCellspaceInfos() const;
      bool delSubCellspaceInfo_glbID(int glbID);
      bool delSubCellspaceInfo_lclID(int lclID);
      void clearSubCellspaceInfos();
      const pRPL::IntVect allSubCellspaceIDs() const;
      int maxSubCellspaceID() const;
      bool calcAllBRs(bool onlyUpdtCtrCell,
                      const pRPL::Neighborhood *pNbrhd);
      
      pRPL::SubCellspaceInfo* subCellspaceInfo_glbID(int glbID,
                                                     bool warning = false);
      const pRPL::SubCellspaceInfo* subCellspaceInfo_glbID(int glbID,
                                                           bool warning = false) const;
      pRPL::SubCellspaceInfo* subCellspaceInfo_lclID(int lclID,
                                                     bool warning = false);
      const pRPL::SubCellspaceInfo* subCellspaceInfo_lclID(int lclID,
                                                           bool warning = false) const;
      
      /* SubCellspace */
      int nSubCellspaces() const;
      pRPL::SubCellspace* addSubCellspace(int glbID,
                                          double *pInitVal = NULL);
      bool delSubCellspace_glbID(int glbID);
      bool delSubCellspace_lclID(int lclID);
      void clearSubCellspaces();
      
      pRPL::SubCellspace* subCellspace_glbID(int glbID,
                                             bool warning = false);
      const pRPL::SubCellspace* subCellspace_glbID(int glbID,
                                                   bool warning = false) const;
      pRPL::SubCellspace* subCellspace_lclID(int lclID,
                                             bool warning = false);
      const pRPL::SubCellspace* subCellspace_lclID(int lclID,
                                                   bool warning = false) const;
      
      void setUpdateTracking(bool toTrack);
      void allUpdatedIdxs(vector<long> &vUpdtGlbIdxs) const;
      void clearUpdateTracks();
      bool loadCellStream(void *aStream, int nCells);

      /*ByGDAL------------*/
      bool loadSubCellspaceByGDAL(int glbID);
      bool writeSubCellspaceByGDAL(int glbID);

      bool writeTmpSubCellspaceByGDAL(int glbID,
                                      int alias = 0);
      bool loadTmpSubCellspaceByGDAL(int glbID,
                                     int alias = 0);
      bool mergeSubCellspaceByGDAL(int glbID,
                                   int alias = 0);
      

      /*-------byPGTIOL--------------------*/
      bool loadSubCellspaceByPGTIOL(int glbID);
      bool writeSubCellspaceByPGTIOL(int glbID);

      /*
      bool writeTmpSubCellspaceByPGTIOL(int glbID,
                                        int alias = 0);
      bool loadTmpSubCellspaceByPGTIOL(int glbID,
                                       int alias = 0);
      bool mergeSubCellspaceByPGTIOL(int glbID,
                                     int alias = 0);
      */

    private:
      const string _tmpSubCellspaceFileName(int glbID,
                                            int alias = 0) const;
      bool _parseDsName(const char *aFileName);

    private:
      string _name;

      pRPL::CellspaceInfo *_pCellspcInfo;
      pRPL::Cellspace *_pCellspc;

      list<pRPL::SubCellspaceInfo> _lSubInfos;
      list<pRPL::CellspaceInfo> _lCellspcInfos;
      list<pRPL::SubCellspace> _lSubCellspcs;

      GDALDataset *_pGdalDS;
	    PGTIOLDataset *_pPgtiolDS;
      int _iBand;
      string _dsPath;
      string _dsName;
  };
};

#endif
