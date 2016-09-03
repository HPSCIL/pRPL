#include "prpl-layer.h"

const string pRPL::Layer::
_tmpSubCellspaceFileName(int glbID,
                         int alias) const {
  string tmpFileName;
  tmpFileName.append("tmp_");
  if(alias > 0) {
    stringstream ssAlias;
    ssAlias << alias;
    tmpFileName.append(ssAlias.str() + "_");
  }
  tmpFileName.append(_name);
  stringstream ssID;
  ssID << glbID;
  tmpFileName.append(ssID.str());
  //tmpFileName.append(".tif");
  return tmpFileName;
}

bool pRPL::Layer::
_parseDsName(const char *aFileName) {

  return true;
}

pRPL::Layer::
Layer(const char *aName)
  :_pCellspcInfo(NULL),
   _pCellspc(NULL),
   _pGdalDS(NULL),
   _pPgtiolDS(NULL),
   _iBand(0) {  
if(aName != NULL) {
   _name.assign(aName);
  }
  else {
    _name.assign("Untitled_Layer");
  }
}

pRPL::Layer::
~Layer() {
  delCellspaceInfo();
  delCellspace();
  clearSubCellspaceInfos();
  clearSubCellspaces();
  closeGdalDS();
  closePgtiolDS();
}

const char* pRPL::Layer::
name() const {
  return _name.c_str();
}

void pRPL::Layer::
name(const char *aName) {
  _name.assign(aName);
}

const pRPL::SpaceDims* pRPL::Layer::
glbDims() const {
  if(_pCellspcInfo == NULL) {
    return NULL;
  }
  return &(_pCellspcInfo->dims());
}

const char* pRPL::Layer::
dataType() const {
  if(_pCellspcInfo == NULL) {
    return NULL;
  }
  return _pCellspcInfo->dataType();
}

size_t pRPL::Layer::
dataSize() const {
  if(_pCellspcInfo == NULL) {
    return 0;
  }
  return _pCellspcInfo->dataSize();
}

long pRPL::Layer::tileSize() const{
  if(_pCellspcInfo == NULL) {
    return 0;
  }
  return _pCellspcInfo->tileSize();
}

const pRPL::CellspaceGeoinfo* pRPL::Layer::
glbGeoinfo() const {
  if(_pCellspcInfo == NULL) {
    return NULL;
  }
  return _pCellspcInfo->georeference();
}

bool pRPL::Layer::
openGdalDS(const char *aFileName,
           GDALAccess ioOption) {
  closeGdalDS();
  if(aFileName == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to open a GDAL dataset using a NULL file name" \
         << endl;
    return false;
  }
  _pGdalDS = (GDALDataset *)GDALOpen(aFileName, ioOption);
  if(_pGdalDS == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to open GDAL file (" \
         << aFileName << ")" << endl;
    return false;
  }
  return true;
}

int pRPL::Layer::
dsBand() const {
  return _iBand;
}

void pRPL::Layer::
dsBand(int iBand) {
  _iBand = iBand;
}

void pRPL::Layer::
closeGdalDS() {
  if(_pGdalDS != NULL) {
    GDALClose((GDALDatasetH)_pGdalDS);
    _pGdalDS = NULL;
  }
  _iBand = 0;
}

bool pRPL::Layer::
createGdalDS(const char *aFileName,
             const char *aGdalFormat,
             char **aGdalOptions) {
  closeGdalDS();

  if(_pCellspcInfo == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: NULL CellspaceInfo, unable to create a GDAL dataset" \
         << endl;
    return false;
  }
  if(_pCellspcInfo->isEmpty(true)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to create a GDAL dataset with an empty CellspaceInfo" \
         << endl;
    return false;
  }
  if(aFileName == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to create a GDAL dataset using a NULL file name" \
         << endl;
    return false;
  }
  if(aGdalFormat == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to create a GDAL dataset using a NULL GDAL format" \
         << endl;
    return false;
  }

  GDALDriver *pGdalDriver = GetGDALDriverManager()->GetDriverByName(aGdalFormat);
  if(pGdalDriver == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to load the GDAL Driver for (" \
         << aGdalFormat << ")" << endl;
    return false;
  }

  GDALDataType gdalDataType;
  pRPL::cppType2gdalType(gdalDataType, string(_pCellspcInfo->dataType()));
  if(gdalDataType == GDT_Unknown) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to determine appropriate data type for the GDAL dataset" \
         << endl;
    return false;
  }

  _pGdalDS = pGdalDriver->Create(aFileName,
                                 _pCellspcInfo->nCols(), _pCellspcInfo->nRows(), 1,
                                 gdalDataType, aGdalOptions);
  if(_pGdalDS == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to create GDAL dataset (" << aFileName << ")" << endl;
    return false;
  }

  _iBand = 1;
  
  _pGdalDS->GetRasterBand(_iBand)->SetNoDataValue(_pCellspcInfo->getNoDataValAs<double>());
  //_pGdalDS->GetRasterBand(_iBand)->Fill(_pCellspcInfo->getNoDataValAs<double>());

  if(_pCellspcInfo->isGeoreferenced(false)) {
    _pGdalDS->SetProjection(_pCellspcInfo->georeference()->projection().c_str());
    double aGeoTransform[6];
    _pCellspcInfo->georeference()->geoTransform(aGeoTransform);
    _pGdalDS->SetGeoTransform(aGeoTransform);
  }

  return true;
}

bool pRPL::Layer::
openPgtiolDS(const char *aFileName,
             int prc,
             PGTIOLAccess ioOption) {
  closePgtiolDS();
  if(aFileName == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to open a GDAL dataset using a NULL file name" \
         << endl;
    return false;
  }
  _pPgtiolDS = new PGTIOLDataset(aFileName, ioOption, prc);

  if(_pPgtiolDS == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to open GDAL file (" \
         << aFileName << ")" << endl;
    return false;
  }
  if(ioOption == PG_Update){
   _pPgtiolDS->SetWriteProcess(prc);
   _pPgtiolDS->SetOutputGTiff();
  }
  return true;
}

void pRPL::Layer::
closePgtiolDS() {
  if(_pPgtiolDS != NULL) {
    delete _pPgtiolDS;
    _pPgtiolDS = NULL;
  }
  _iBand = 0;
}

bool pRPL::Layer::
createPgtiolDS(const char *aFileName,
               int prc,
               char **aPgtOptions) {
  
  closePgtiolDS();

  if(_pCellspcInfo == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: NULL CellspaceInfo, unable to create a GDAL dataset" \
         << endl;
    return false;
  }
  if(_pCellspcInfo->isEmpty(true)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to create a GDAL dataset with an empty CellspaceInfo" \
         << endl;
    return false;
  }
  if(aFileName == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to create a GDAL dataset using a NULL file name" \
         << endl;
    return false;
  }

  GDALDataType gdalDataType;
  pRPL::cppType2gdalType(gdalDataType, string(_pCellspcInfo->dataType()));
  if(gdalDataType == GDT_Unknown) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to determine appropriate data type for the GDAL dataset" \
         << endl;
    return false;
  }

  _pPgtiolDS = new PGTIOLDataset(aFileName, PG_Update, aPgtOptions,
                                 _pCellspcInfo->nCols(),_pCellspcInfo->nRows(), 1,gdalDataType);

  if(_pPgtiolDS == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to create pGTIOL dataset (" << aFileName << ")" << endl;
    return false;
  }

  _iBand = 1;
  _pPgtiolDS->SetRasterBand(_iBand);
  _pPgtiolDS->SetNoDataValue(_pCellspcInfo->getNoDataValAs<double>());
  if(_pCellspcInfo->isGeoreferenced(false)) {
	  _pPgtiolDS->SetProjection(_pCellspcInfo->georeference()->projection());
     double aGeoTransform[6];
     _pCellspcInfo->georeference()->geoTransform(aGeoTransform);
     _pPgtiolDS->SetGeoTransform(aGeoTransform);
	 
	 _pPgtiolDS->SetTileWidth(_pCellspcInfo->tileSize());
	 _pPgtiolDS->SetTileLength(_pCellspcInfo->tileSize());
  }
   _pPgtiolDS->SetWriteProcess(prc);
   _pPgtiolDS->SetOutputGTiff();
   return true;
}

bool pRPL::Layer::
decompose(int nRowSubspcs,
          int nColSubspcs,
          const pRPL::Neighborhood *pNbrhd) {
  clearSubCellspaceInfos();
  clearSubCellspaces();

  if(_pCellspcInfo == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: CellspaceInfo has NOT been initialized to be decomposed" \
         << endl;
    return false;
  }

  pRPL::DomDcmpMethod dcmpMethod = pRPL::SNGL_PRC;
  if(nRowSubspcs == 1 && nColSubspcs == 1) {
    dcmpMethod = pRPL::SNGL_PRC;
  }
  else if(nRowSubspcs > 1 && nColSubspcs == 1) {
    dcmpMethod = pRPL::SMPL_ROW;
  }
  else if(nRowSubspcs == 1 && nColSubspcs > 1) {
    dcmpMethod = pRPL::SMPL_COL;
  }
  else if(nRowSubspcs > 1 && nColSubspcs > 1) {
    dcmpMethod = pRPL::SMPL_BLK;
  }
  else {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid numbers of SubCellspaces (" \
         << nRowSubspcs << ", " << nColSubspcs << ")" \
         << endl;
    return false;
  }

  pRPL::SmplDcmp dcmp;
  switch(dcmpMethod) {
    case pRPL::SNGL_PRC:
      if(!dcmp.rowDcmp(_lSubInfos, *glbDims(), nRowSubspcs, pNbrhd)) {
        return false;
      }
      break;
    case pRPL::SMPL_ROW:
      if(!dcmp.rowDcmp(_lSubInfos, *glbDims(), nRowSubspcs, pNbrhd)) {
        return false;
      }
      break;
    case pRPL::SMPL_COL:
      if(!dcmp.colDcmp(_lSubInfos, *glbDims(), nColSubspcs, pNbrhd)) {
        return false;
      }
      break;
    case pRPL::SMPL_BLK:
      if(!dcmp.blkDcmp(_lSubInfos, *glbDims(), nRowSubspcs, nColSubspcs, pNbrhd)) {
        return false;
      }
      break;
    default:
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: the decomposition method (" << dcmpMethod \
           << ") is not supported." << endl;
      break;
  }

  return true;
}

bool pRPL::Layer::
copyDcmp(const Layer *pFromLyr) {
  clearSubCellspaceInfos();
  clearSubCellspaces();

  if(*glbDims() != *(pFromLyr->glbDims())) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unequal global spatial dimensions: (" \
         << *glbDims() << ") and (" << *(pFromLyr->glbDims()) \
         << ")" << endl;
    return false;
  }

  for(int iSubInfo = 0; iSubInfo < pFromLyr->nSubCellspaceInfos(); iSubInfo++) {
    _lSubInfos.push_back(pRPL::SubCellspaceInfo(*(pFromLyr->subCellspaceInfo_lclID(iSubInfo))));
  }
  return true;
}

bool pRPL::Layer::
isDecomposed() const {
  return (!_lSubInfos.empty());
}

void pRPL::Layer::
initCellspaceInfo() {
  delCellspaceInfo();
  delCellspace();
  _pCellspcInfo = new pRPL::CellspaceInfo();
}

void pRPL::Layer::
initCellspaceInfo(const pRPL::SpaceDims &dims,
                  const char *aTypeName,
                  size_t typeSize,
                  const pRPL::CellspaceGeoinfo *pGeoinfo,
                  long tileSize) {
  delCellspaceInfo();
  delCellspace();
  _pCellspcInfo = new pRPL::CellspaceInfo(dims, aTypeName, typeSize, pGeoinfo, tileSize);
}

void pRPL::Layer::
initCellspaceInfo(const pRPL::CellspaceInfo &cellspcInfo) {
  delCellspaceInfo();
  delCellspace();
  _pCellspcInfo = new pRPL::CellspaceInfo(cellspcInfo);
}

/*----------------GDAL-----------------*/
bool pRPL::Layer::
initCellspaceInfoByGDAL() {
  delCellspaceInfo();
  delCellspace();

  _pCellspcInfo = new pRPL::CellspaceInfo();
  if(!_pCellspcInfo->initByGDAL(_pGdalDS, _iBand, true)) {
    delCellspaceInfo();
    return false;
  }
  return true;
}

/*----------------PGTIOL-----------------*/
bool pRPL::Layer::
initCellspaceInfoByPGTIOL() {
  delCellspaceInfo();
  delCellspace();

  _pCellspcInfo = new pRPL::CellspaceInfo();
  if(!_pCellspcInfo->initByPGTIOL(_pPgtiolDS, _iBand, true)) {
    delCellspaceInfo();
    return false;
  }
  return true;
}

void pRPL::Layer::
delCellspaceInfo() {
   if(_pCellspcInfo != NULL) {
    delete _pCellspcInfo;
    _pCellspcInfo = NULL;
  }
}

pRPL::CellspaceInfo* pRPL::Layer::
cellspaceInfo() {
  return _pCellspcInfo;
}

const pRPL::CellspaceInfo* pRPL::Layer::
cellspaceInfo() const {
  return _pCellspcInfo;
}

bool pRPL::Layer::
createCellspace(double *pInitVal) {
  delCellspace();
  if(_pCellspcInfo == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: NULL CellspaceInfo to create a Cellspace" \
         << endl;
    return false;
  }
  _pCellspc = new pRPL::Cellspace(_pCellspcInfo, true);

  if(pInitVal != NULL) {
    if(!_pCellspc->initValAs<double>(pInitVal, true)) {
      return false;
    }
  }

  return true;
}


/*----------------GDAL--------------------*/
bool pRPL::Layer::
loadCellspaceByGDAL() {
  if(!initCellspaceInfoByGDAL()) {
    return false;
  }
  if(!createCellspace()) {
    return false;
  }
  return _pCellspc->readGDALData(_pGdalDS, _iBand, NULL, true);
}

bool pRPL::Layer::
writeCellspaceByGDAL() {
  if(_pCellspcInfo == NULL ||
     _pCellspc == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: NULL CellspaceInfo or NULL Cellspace, unable to write to a GDAL dataset" \
         << endl;
    return false;
  }

  return _pCellspc->writeGDALData(_pGdalDS, _iBand, NULL, NULL, true);
}

/*---------------------PGTIOL----------------------------*/
bool pRPL::Layer::loadCellspaceByPGTIOL() {
 if(!initCellspaceInfoByPGTIOL()) {
	return false;
  }
  if(!createCellspace()) {
	return false;
  }
  return _pCellspc->readPGTIOLData(_pPgtiolDS, _iBand, NULL, true);
}

bool pRPL::Layer::writeCellspaceByPGTIOL() {
  if(_pCellspcInfo == NULL ||
	 _pCellspc == NULL) {
	cerr << __FILE__ << " " << __FUNCTION__ \
		 << " Error: NULL CellspaceInfo or NULL Cellspace, unable to write to a GDAL dataset" \
		 << endl;
	return false;
  }
  
  return _pCellspc->writePGTIOLData(_pPgtiolDS, _iBand, NULL, NULL, true);
}

void pRPL::Layer::
delCellspace() {
  if(_pCellspc) {
    delete _pCellspc;
    _pCellspc = NULL;
  }
}

pRPL::Cellspace* pRPL::Layer::
cellspace() {
  return _pCellspc;
}

const pRPL::Cellspace* pRPL::Layer::
cellspace() const {
  return _pCellspc;
}

int pRPL::Layer::
nSubCellspaceInfos() const {
  return _lSubInfos.size();
}

const list<pRPL::SubCellspaceInfo>* pRPL::Layer::
allSubCellspaceInfos() const {
  return &_lSubInfos;
}

bool pRPL::Layer::
delSubCellspaceInfo_glbID(int glbID) {
  bool found = false;
  list<pRPL::SubCellspaceInfo>::iterator itrSubInfo = _lSubInfos.begin();
  while(itrSubInfo != _lSubInfos.end()) {
    if(glbID == itrSubInfo->id()) {
      _lSubInfos.erase(itrSubInfo);
      found = true;
      break;
    }
    itrSubInfo++;
  }
  return found;
}

bool pRPL::Layer::
delSubCellspaceInfo_lclID(int lclID) {
  bool found = false;
  if(lclID >= 0 && lclID < _lSubInfos.size()) {
    list<pRPL::SubCellspaceInfo>::iterator itrSubInfo = _lSubInfos.begin();
    advance(itrSubInfo, lclID);
    _lSubInfos.erase(itrSubInfo);
    found = true;
  }
  return found;
}
      
void pRPL::Layer::
clearSubCellspaceInfos() {
  _lSubInfos.clear();
}

pRPL::SubCellspaceInfo* pRPL::Layer::
subCellspaceInfo_glbID(int glbID,
                       bool warning) {
  pRPL::SubCellspaceInfo *pSubInfo = NULL;
  list<pRPL::SubCellspaceInfo>::iterator itrSubInfo = _lSubInfos.begin();
  while(itrSubInfo != _lSubInfos.end()) {
    if(glbID == itrSubInfo->id()) {
      pSubInfo = &(*itrSubInfo);
      break;
    }
    itrSubInfo++;
  }
  if(pSubInfo == NULL && warning) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Warning: unable to find the SubCellspaceInfo with the global ID (" \
         << glbID << ")" << endl;
  }
  return pSubInfo;
}

const pRPL::SubCellspaceInfo* pRPL::Layer::
subCellspaceInfo_glbID(int glbID,
                       bool warning) const {
  const pRPL::SubCellspaceInfo *pSubInfo = NULL;
  list<pRPL::SubCellspaceInfo>::const_iterator itrSubInfo = _lSubInfos.begin();
  while(itrSubInfo != _lSubInfos.end()) {
    if(glbID == itrSubInfo->id()) {
      pSubInfo = &(*itrSubInfo);
      break;
    }
    itrSubInfo++;
  }
  if(pSubInfo == NULL && warning) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Warning: unable to find the SubCellspaceInfo with the global ID (" \
        << glbID << ")" << endl;
  }
  return pSubInfo;
}

pRPL::SubCellspaceInfo* pRPL::Layer::
subCellspaceInfo_lclID(int lclID,
                       bool warning) {
  pRPL::SubCellspaceInfo *pSubInfo = NULL;
  if(lclID >= 0 && lclID < _lSubInfos.size()) {
    list<pRPL::SubCellspaceInfo>::iterator itrSubInfo = _lSubInfos.begin();
    advance(itrSubInfo, lclID);
    pSubInfo = &(*itrSubInfo);
  }
  if(pSubInfo == NULL && warning) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Warning: unable to find the SubCellspaceInfo with the local ID (" \
        << lclID << ")" << endl;
  }
  return pSubInfo;
}

const pRPL::SubCellspaceInfo* pRPL::Layer::
subCellspaceInfo_lclID(int lclID,
                       bool warning) const {
  const pRPL::SubCellspaceInfo *pSubInfo = NULL;
  if(lclID >= 0 && lclID < _lSubInfos.size()) {
    list<pRPL::SubCellspaceInfo>::const_iterator itrSubInfo = _lSubInfos.begin();
    advance(itrSubInfo, lclID);
    pSubInfo = &(*itrSubInfo);
  }
  if(pSubInfo == NULL && warning) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Warning: unable to find the SubCellspaceInfo with the local ID (" \
        << lclID << ")" << endl;
  }
  return pSubInfo;
}

const pRPL::IntVect pRPL::Layer::
allSubCellspaceIDs() const {
  pRPL::IntVect vSubspcIDs;
  list<pRPL::SubCellspaceInfo>::const_iterator itrSubInfo = _lSubInfos.begin();
  while(itrSubInfo != _lSubInfos.end()) {
    vSubspcIDs.push_back(itrSubInfo->id());
    itrSubInfo++;
  }
  return vSubspcIDs;
}

int pRPL::Layer::
maxSubCellspaceID() const {
  int maxSubID = pRPL::ERROR_ID;
  list<pRPL::SubCellspaceInfo>::const_iterator itrSubInfo = _lSubInfos.begin();
  while(itrSubInfo != _lSubInfos.end()) {
    if(maxSubID < itrSubInfo->id()) {
      maxSubID = itrSubInfo->id();
    }
    itrSubInfo++;
  }
  return maxSubID;
}

bool pRPL::Layer::
calcAllBRs(bool onlyUpdtCtrCell,
           const pRPL::Neighborhood *pNbrhd) {
  list<pRPL::SubCellspaceInfo>::iterator itrSubInfo = _lSubInfos.begin();
  while(itrSubInfo != _lSubInfos.end()) {
    if(!itrSubInfo->calcBRs(onlyUpdtCtrCell, pNbrhd, &_lSubInfos)) {
      return false;
    }
    itrSubInfo++;
  }
  return true;
}

int pRPL::Layer::
nSubCellspaces() const {
  return _lSubCellspcs.size();
}

pRPL::SubCellspace* pRPL::Layer::
addSubCellspace(int glbID,
                double *pInitVal) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID);
  if(pSubspc != NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubCellspace with the same global ID (" \
         << glbID << ") already exists. Failed to add SubCellspace." \
         << endl;
  }
  else {
    pRPL::SubCellspaceInfo* pSubspcInfo = subCellspaceInfo_glbID(glbID);
    if(pSubspcInfo == NULL) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error: unable to find SubCellspaceInfo with ID (" \
          << glbID << ")" << endl;
    }
    else if(_pCellspcInfo == NULL) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error: CellspaceInfo in this Layer has NOT been initialized yet" \
          << endl;
    }
    else if(std::strcmp(dataType(), "UNKNOWN") == 0 || dataSize() <= 0) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error: invalid data type (" << dataType() \
          << ") or data size (" << dataSize() << ")" << endl;
    }
    else {
      _lCellspcInfos.push_back(pRPL::CellspaceInfo(pSubspcInfo->dims(),dataType(),dataSize()));
      pRPL::CellspaceInfo *pCellspcInfo_sub = &(_lCellspcInfos.back());

      double noDataVal = _pCellspcInfo->getNoDataValAs<double>();
      pCellspcInfo_sub->setNoDataValAs<double>(noDataVal);

      if(_pCellspcInfo->isGeoreferenced(false)) {
        pRPL::CellspaceGeoinfo subspcGeoinfo = _pCellspcInfo->georeference()->subSpaceGeoinfo(*pSubspcInfo);
        pCellspcInfo_sub->georeference(subspcGeoinfo);
      }

      pCellspcInfo_sub->setTileSize(tileSize());

      _lSubCellspcs.push_back(pRPL::SubCellspace(pCellspcInfo_sub, pSubspcInfo));
      pSubspc = &(_lSubCellspcs.back());

      if(pInitVal != NULL) {
        pSubspc->initValAs<double>(pInitVal, true);
        //cout << pSubspc->subInfo()->id() << " initialized noData" << endl;
      }
    }
  }

  return pSubspc;
}

bool pRPL::Layer::
delSubCellspace_glbID(int glbID) {
  bool found = false;
  list<pRPL::SubCellspace>::iterator itrSubspc = _lSubCellspcs.begin();
  list<pRPL::CellspaceInfo>::iterator itrCellspcInfo = _lCellspcInfos.begin();
  while(itrSubspc != _lSubCellspcs.end()) {
    if(glbID == itrSubspc->subInfo()->id()) {
      _lSubCellspcs.erase(itrSubspc);
      _lCellspcInfos.erase(itrCellspcInfo);
      found = true;
      break;
    }
    itrSubspc++;
    itrCellspcInfo++;
  }
  return found;
}

bool pRPL::Layer::
delSubCellspace_lclID(int lclID) {
  bool found = false;
  if(lclID >= 0 && lclID < _lSubCellspcs.size()) {
    list<pRPL::CellspaceInfo>::iterator itrCellspcInfo = _lCellspcInfos.begin();
    list<pRPL::SubCellspace>::iterator itrSubspc = _lSubCellspcs.begin();
    advance(itrCellspcInfo, lclID);
    advance(itrSubspc, lclID);
    _lSubCellspcs.erase(itrSubspc);
    _lCellspcInfos.erase(itrCellspcInfo);
    found = true;
  }
  return found;
}

void pRPL::Layer::
clearSubCellspaces() {
  _lCellspcInfos.clear();
  _lSubCellspcs.clear();
}

pRPL::SubCellspace* pRPL::Layer::
subCellspace_glbID(int glbID,
                   bool warning) {
  pRPL::SubCellspace *pSubspc = NULL;
  list<pRPL::SubCellspace>::iterator itrSubspc = _lSubCellspcs.begin();
  while(itrSubspc != _lSubCellspcs.end()) {
    if(glbID == itrSubspc->subInfo()->id()) {
      pSubspc = &(*itrSubspc);
      break;
    }
    itrSubspc++;
  }
  if(pSubspc == NULL && warning) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Warning: unable to find the SubCellspace with the global ID (" \
         << glbID << ") within Layer (" << _name << ")" << endl;
  }
  return pSubspc;
}

const pRPL::SubCellspace* pRPL::Layer::
subCellspace_glbID(int glbID,
                   bool warning) const {
  const pRPL::SubCellspace *pSubspc = NULL;
  list<pRPL::SubCellspace>::const_iterator itrSubspc = _lSubCellspcs.begin();
  while(itrSubspc != _lSubCellspcs.end()) {
    if(glbID == itrSubspc->subInfo()->id()) {
      pSubspc = &(*itrSubspc);
      break;
    }
    itrSubspc++;
  }
  if(pSubspc == NULL && warning) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Warning: unable to find the SubCellspace with the global ID (" \
         << glbID << ") within Layer (" << _name << ")" << endl;
  }
  return pSubspc;
}

pRPL::SubCellspace* pRPL::Layer::
subCellspace_lclID(int lclID,
                   bool warning) {
  pRPL::SubCellspace *pSubspc = NULL;
  if(lclID >= 0 && lclID < _lSubCellspcs.size()) {
    list<pRPL::SubCellspace>::iterator itrSubspc = _lSubCellspcs.begin();
    advance(itrSubspc, lclID);
    pSubspc = &(*itrSubspc);
  }
  if(pSubspc == NULL && warning) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Warning: unable to find the SubCellspace with the local ID (" \
         << lclID << ")" << endl;
  }
  return pSubspc;
}

const pRPL::SubCellspace* pRPL::Layer::
subCellspace_lclID(int lclID,
                   bool warning) const {
  const pRPL::SubCellspace *pSubspc = NULL;
  if(lclID >= 0 && lclID < _lSubCellspcs.size()) {
    list<pRPL::SubCellspace>::const_iterator itrSubspc = _lSubCellspcs.begin();
    advance(itrSubspc, lclID);
    pSubspc = &(*itrSubspc);
  }
  if(pSubspc == NULL && warning) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Warning: unable to find the SubCellspace with the local ID (" \
         << lclID << ")" << endl;
  }
  return pSubspc;
}

void pRPL::Layer::
setUpdateTracking(bool toTrack) {
  if(_lSubCellspcs.size() > 0) {
    list<pRPL::SubCellspace>::iterator itrSubspc = _lSubCellspcs.begin();
    while(itrSubspc != _lSubCellspcs.end()) {
      itrSubspc->setUpdateTracking(toTrack);
      itrSubspc++;
    }
  }
  else if(_pCellspc !=  NULL && !_pCellspc->isEmpty()) {
    _pCellspc->setUpdateTracking(toTrack);
  }
}

void pRPL::Layer::
allUpdatedIdxs(vector<long> &vUpdtGlbIdxs) const {
  vUpdtGlbIdxs.clear();
  if(_lSubCellspcs.size() > 0) {
    list<pRPL::SubCellspace>::const_iterator itrSubspc = _lSubCellspcs.begin();
    while(itrSubspc != _lSubCellspcs.end()) {
      const pRPL::LongVect &vLclIdxs = itrSubspc->getUpdatedIdxs();
      for(int iIdx = 0; iIdx < vLclIdxs.size(); iIdx++) {
        vUpdtGlbIdxs.push_back(itrSubspc->subInfo()->lclIdx2glbIdx(vLclIdxs[iIdx]));
      }
      itrSubspc++;
    }
  }
  else if(_pCellspc != NULL && !_pCellspc->isEmpty()) {
    vUpdtGlbIdxs = _pCellspc->getUpdatedIdxs();
  }
}

void pRPL::Layer::
clearUpdateTracks() {
  if(_lSubCellspcs.size() > 0) {
    list<pRPL::SubCellspace>::iterator itrSubspc = _lSubCellspcs.begin();
    while(itrSubspc != _lSubCellspcs.end()) {
      itrSubspc->clearUpdatedIdxs();
      itrSubspc++;
    }
  }
  if(_pCellspc != NULL && !_pCellspc->isEmpty()) {
    _pCellspc->clearUpdatedIdxs();
  }
}

bool pRPL::Layer::
loadCellStream(void *aStream,
               int nCells) {
  if(aStream == NULL && nCells > 0) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: NULL stream to load" << endl;
    return false;
  }
  if(nCells < 0) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid number of cells (" << nCells \
         << ")" << endl;
    return false;
  }
  if(_pCellspcInfo == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the Layer has NOT been initialized yet" \
         << endl;
    return false;
  }

  size_t dataSize = _pCellspcInfo->dataSize();
  void *aCellVal = (void *)malloc(dataSize);
  long cellGlbIdx, cellLclIdx;
  void *pCell;

  for(int iCell = 0; iCell < nCells; iCell++) {
    memcpy((char *)(&cellGlbIdx), (char *)aStream+iCell*(sizeof(long) + dataSize), sizeof(long));
    memcpy((char *)aCellVal, (char *)aStream+iCell*(sizeof(long) + dataSize)+sizeof(long), dataSize);
    list<pRPL::SubCellspace>::iterator itrSubspc = _lSubCellspcs.begin();
    while(itrSubspc != _lSubCellspcs.end()) {
      if(itrSubspc->subInfo()->MBR().ifContain(cellGlbIdx, itrSubspc->subInfo()->glbDims())) {
        cellLclIdx = itrSubspc->subInfo()->glbIdx2lclIdx(cellGlbIdx);
        pCell = itrSubspc->at(cellLclIdx, true);
        if(pCell == NULL) {
          return false;
        }
        memcpy(pCell, aCellVal, dataSize);
        //cout << _name << " [" << cellGlbIdx << ", " << itrSubspc->atAs<unsigned short>(cellLclIdx) << "]" << endl;
      }
      itrSubspc++;
    }
  } // end -- for(iCell)

  free(aCellVal);

  return true;
}

bool pRPL::Layer::
loadSubCellspaceByGDAL(int glbID) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID, false);
  if(pSubspc == NULL) {
    pSubspc = addSubCellspace(glbID);
    if(pSubspc == NULL) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: unable to create a SubCellspace with the global ID (" \
           << glbID << ")" << endl;
      return false;
    }
  }

  const CellCoord gdalOff = pSubspc->subInfo()->MBR().nwCorner();
  return pSubspc->readGDALData(_pGdalDS, _iBand, &gdalOff, true);
}

bool pRPL::Layer::
writeSubCellspaceByGDAL(int glbID) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID, true);
  if(pSubspc == NULL) {
    return false;
  }

  const pRPL::CellCoord gdalOff = pSubspc->subInfo()->MBR().nwCorner();
  const pRPL::CoordBR workBR = pSubspc->subInfo()->workBR();
  return pSubspc->writeGDALData(_pGdalDS, _iBand, &gdalOff, &workBR, true);
}

bool pRPL::Layer::
writeTmpSubCellspaceByGDAL(int glbID,
                           int alias) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID, true);
  if(pSubspc == NULL) {
    return false;
  }

  if(_pGdalDS == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the layer's GDAL dataset has not been initialized" \
         << endl;
    return false;
  }
  GDALDriver *pGdalDriver = _pGdalDS->GetDriver();

  GDALDataset *pTmpGdalDS;
  string tmpFileName = _tmpSubCellspaceFileName(glbID, alias);
  if(!pSubspc->createGdalDS(pTmpGdalDS, tmpFileName.c_str(), pGdalDriver)) {
    return false;
  }

  if(!pSubspc->writeGDALData(pTmpGdalDS)) {
    return false;
  }

  GDALClose((GDALDatasetH)pTmpGdalDS);

  return true;
}

bool pRPL::Layer::
loadTmpSubCellspaceByGDAL(int glbID,
                          int alias) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID, false);
  if(pSubspc == NULL) {
    pSubspc = addSubCellspace(glbID);
    if(pSubspc == NULL) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error: unable to create a SubCellspace with the global ID (" \
          << glbID << ")" << endl;
      return false;
    }
  }

  string tmpFileName = _tmpSubCellspaceFileName(glbID, alias);
  GDALDataset *pTmpGdalDS;
  pTmpGdalDS = (GDALDataset *)GDALOpen(tmpFileName.c_str(), GA_ReadOnly);
  if(pTmpGdalDS == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to open GDAL file (" \
         << tmpFileName << ")" << endl;
    return false;
  }

  if(!pSubspc->readGDALData(pTmpGdalDS)) {
    return false;
  }

  GDALClose((GDALDatasetH)pTmpGdalDS);

  // delete the temporary file
  GDALDriver *pGdalDriver = _pGdalDS->GetDriver();
  if(pGdalDriver->Delete(_tmpSubCellspaceFileName(glbID, alias).c_str()) == CE_Failure) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: failed to delete temporary dataset ( " \
        << _tmpSubCellspaceFileName(glbID, alias) << ")" \
        << endl;
    return false;
  }

  return true;
}

bool pRPL::Layer::
mergeSubCellspaceByGDAL(int glbID,
                        int alias) {
  if(!loadTmpSubCellspaceByGDAL(glbID, alias) ||
     !writeSubCellspaceByGDAL(glbID)) {
    return false;
  }

  return true;
}


/*-------byPGTIOL--------------------*/
bool pRPL::Layer::
loadSubCellspaceByPGTIOL(int glbID) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID, false);
  if(pSubspc == NULL) {
    pSubspc = addSubCellspace(glbID);
    if(pSubspc == NULL) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: unable to create a SubCellspace with the global ID (" \
           << glbID << ")" << endl;
      return false;
    }
  }

  const CellCoord gdalOff = pSubspc->subInfo()->MBR().nwCorner();
  return pSubspc->readPGTIOLData(_pPgtiolDS, _iBand, &gdalOff, true);
}

bool pRPL::Layer::
writeSubCellspaceByPGTIOL(int glbID) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID, true);
  if(pSubspc == NULL) {
    return false;
  }
  
  const pRPL::CellCoord gdalOff = pSubspc->subInfo()->MBR().nwCorner();
  const pRPL::CoordBR workBR = pSubspc->subInfo()->workBR();
  return pSubspc->writePGTIOLData(_pPgtiolDS, _iBand, &gdalOff, &workBR, true);
}

/*
bool pRPL::Layer::
writeTmpSubCellspaceByPGTIOL(int glbID,
                           int alias) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID, true);
  if(pSubspc == NULL) {
    return false;
  }

  if(_pGdalDS == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the layer's GDAL dataset has not been initialized" \
         << endl;
    return false;
  }
  PGTIOLDriver *pPgtiolDriver = _pPgtiolDS->GetDriver();

  PGTIOLDataset *pTmpPgtiolDS;
  string tmpFileName = _tmpSubCellspaceFileName(glbID, alias);
  if(!pSubspc->createPgtiolDS(pTmpPgtiolDS, tmpFileName.c_str(), pPgtiolDriver)) {
    return false;
  }

  if(!pSubspc->writePGTIOLData(pTmpPgtiolDS)) {
    return false;
  }

  PGTIOLClose((PGTIOLDataset*)pTmpPgtiolDS);

  return true;
}

bool pRPL::Layer::
loadTmpSubCellspaceByPGTIOL(int glbID,
                          int alias) {
  pRPL::SubCellspace *pSubspc = subCellspace_glbID(glbID, false);
  if(pSubspc == NULL) {
    pSubspc = addSubCellspace(glbID);
    if(pSubspc == NULL) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error: unable to create a SubCellspace with the global ID (" \
          << glbID << ")" << endl;
      return false;
    }
  }

  string tmpFileName = _tmpSubCellspaceFileName(glbID, alias);
  PGTIOLDataset *pTmpPgtiolDS;
  pTmpPgtiolDS = (PGTIOLDataset *)PGTIOLOpen(const_cast<char*>(tmpFileName.c_str()), PG_ReadOnly);
  if(pTmpPgtiolDS == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to open GDAL file (" \
         << tmpFileName << ")" << endl;
    return false;
  }

  if(!pSubspc->readPGTIOLData(pTmpPgtiolDS)) {
    return false;
  }

  PGTIOLClose((PGTIOLDataset*)pTmpPgtiolDS);

  // delete the temporary file
  PGTIOLDriver *pPgtiolDriver = _pPgtiolDS->GetDriver();
  if(pPgtiolDriver->Delete(_tmpSubCellspaceFileName(glbID, alias).c_str()) == PG_Failure) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: failed to delete temporary dataset ( " \
        << _tmpSubCellspaceFileName(glbID, alias) << ")" \
        << endl;
    return false;
  }

  return true;
}

bool pRPL::Layer::
mergeSubCellspaceByPGTIOL(int glbID,
                        int alias) {
  if(!loadTmpSubCellspaceByPGTIOL(glbID, alias) ||
     !writeSubCellspaceByPGTIOL(glbID)) {
    return false;
  }

  return true;
}
*/
