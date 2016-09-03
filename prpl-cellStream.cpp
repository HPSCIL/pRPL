#include "prpl-cellStream.h"

vector<pair<unsigned int, unsigned int> >::const_iterator pRPL::CellStream::
_getLayerInfo(int lyrID) const {
  vector<pair<unsigned int, unsigned int> >::const_iterator itrInfo = _vLayerInfos.begin();
  while(itrInfo != _vLayerInfos.end()) {
    if(itrInfo->first == lyrID) {
      break;
    }
    itrInfo++;
  }
  return itrInfo;
}

pRPL::CellStream::
CellStream()
  :_size(0),
   _aStream(NULL) {}

pRPL::CellStream::
~CellStream() {
  clear();
}

pRPL::CellStream& pRPL::CellStream::
operator=(const pRPL::CellStream &rhs) {
  if(this != &rhs) {
    clear();
    _size = rhs._size;
    _vCellCounts = rhs._vCellCounts;
    _vLayerInfos = rhs._vLayerInfos;
    _aStream = (void *)malloc(_size);
    memcpy((char *)_aStream, (char *)(rhs._aStream), _size);
  }
  return *this;
}

void pRPL::CellStream::
clear() {
  if(_aStream != NULL) {
    free(_aStream);
    _aStream = NULL;
  }
  _size = 0;
  _vCellCounts.clear();
  _vLayerInfos.clear();
}

size_t pRPL::CellStream::
size() const {
  return _size;
}

int pRPL::CellStream::
getTotalCellCount() const {
  int nCells = 0;
  for(int iLyr = 0; iLyr < _vCellCounts.size(); iLyr++) {
    nCells += _vCellCounts[iLyr];
  }
  return nCells;
}

vector<int>& pRPL::CellStream::
getCellCounts() {
  return _vCellCounts;
}

const vector<pair<unsigned int, unsigned int> >& pRPL::CellStream::
getLayerInfos() const {
  return _vLayerInfos;
}

bool pRPL::CellStream::
resize() {
  if(_aStream != NULL) {
    free(_aStream);
    _aStream = NULL;
  }
  _size = 0;

  for(int iLyr = 0; iLyr < _vLayerInfos.size(); iLyr++) {
    unsigned int dataSize = _vLayerInfos[iLyr].second;
    int nCells = _vCellCounts[iLyr];
    if(dataSize <= 0 || nCells < 0) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: invalid data size (" << dataSize \
           << "), or invalid number of cells (" << nCells \
           << "), for Layer (" << _vLayerInfos[iLyr].first \
           << ")" << endl;
      return false;
    }
    _size += (sizeof(long) + dataSize) * nCells;
  }

  _aStream = (void *)malloc(_size);

  return true;
}

const vector<int>& pRPL::CellStream::
getCellCounts() const {
  return _vCellCounts;
}

void* pRPL::CellStream::
getStream() {
  return _aStream;
}

const void* pRPL::CellStream::
getStream() const {
  return _aStream;
}

bool pRPL::CellStream::
addLayer(int lyrID,
         unsigned int typeSize) {
  if(lyrID < 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid Layer ID (" << lyrID \
         << ")" << endl;
    return false;
  }

  if(typeSize <= 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid data size (" << typeSize \
         << ")" << endl;
    return false;
  }

  _vLayerInfos.push_back(make_pair(lyrID, typeSize));
  _vCellCounts.push_back(0);

  return true;
}

bool pRPL::CellStream::
addCell(long cellIdx,
        const void *aCellVal) {
  if(cellIdx < 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: invalid cell index (" << cellIdx << ")" \
        << endl;
    return false;
  }
  if(aCellVal == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: NULL pointer to the cell's value" \
        << endl;
    return false;
  }

  unsigned int typeSize = _vLayerInfos.back().second;
  _aStream = (char *)realloc(_aStream, _size+sizeof(long)+typeSize);
  if(!_aStream) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: failed to re-allocate memory of the CellStream" << endl;
    return false;
  }

  memcpy((char *)((char *)_aStream+_size), &cellIdx, sizeof(long));
  memcpy((char *)((char *)_aStream+_size+sizeof(long)), aCellVal, typeSize);

  _size += sizeof(long) + typeSize;
  _vCellCounts.back()++;

  return true;
}

unsigned int pRPL::CellStream::
getNumLayers() const {
  return _vCellCounts.size();
}

unsigned int pRPL::CellStream::
getTypeSizeOnLayer(int lyrID) const {
  vector<pair<unsigned int, unsigned int> >::const_iterator itrInfo = _getLayerInfo(lyrID);
  return (itrInfo == _vLayerInfos.end()) ? 0 : itrInfo->second;
}

unsigned int pRPL::CellStream::
getNumCellsOnLayer(int lyrID) const {
  vector<pair<unsigned int, unsigned int> >::const_iterator itrInfo = _getLayerInfo(lyrID);
  return (itrInfo == _vLayerInfos.end()) ? 0 : _vCellCounts[itrInfo - _vLayerInfos.begin()];
}

void* pRPL::CellStream::
getCellsOnLayer(int lyrID) {
  void *pCells = NULL;

  vector<pair<unsigned int, unsigned int> >::const_iterator itrInfo = _getLayerInfo(lyrID);
  if(itrInfo != _vLayerInfos.end()) {
    size_t nBytesBefore = 0;
    int iItem = 0;
    unsigned int iLyr2, typeSize2, nCells2;
    vector<pair<unsigned int, unsigned int> >::const_iterator itrInfo2 = _vLayerInfos.begin();
    while(itrInfo2 < itrInfo) {
      typeSize2 = itrInfo2->second;
      nCells2 = _vCellCounts[itrInfo2 - _vLayerInfos.begin()];
      nBytesBefore += (sizeof(long) + typeSize2) * nCells2;
      itrInfo2++;
    } // end -- while(itrInfo2 < itrInfo)

    pCells = (char *)_aStream + nBytesBefore;
  } // end -- if(itrInfo != _vLayerInfos.end())

  return pCells;
}
