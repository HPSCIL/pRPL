#include "prpl-subCellspace.h"

pRPL::SubCellspace::
SubCellspace()
  :_pSubInfo(NULL),
   pRPL::Cellspace() {}

pRPL::SubCellspace::
SubCellspace(pRPL::CellspaceInfo *pInfo,
             pRPL::SubCellspaceInfo *pSubInfo)
  :_pSubInfo(pSubInfo),
   pRPL::Cellspace(pInfo) {}

pRPL::SubCellspace::
SubCellspace(const pRPL::SubCellspace &rhs)
  :_pSubInfo(rhs._pSubInfo),
   pRPL::Cellspace(rhs) {}

pRPL::SubCellspace::
~SubCellspace() {
  _pSubInfo = NULL;
}

pRPL::SubCellspace& pRPL::SubCellspace::
operator=(const pRPL::SubCellspace &rhs) {
  if(this != &rhs) {
    _pSubInfo = rhs._pSubInfo;
    Cellspace::operator=(rhs);
  }
  return *this;
}

pRPL::SubCellspaceInfo* pRPL::SubCellspace::
subInfo() {
  return _pSubInfo;
}

const pRPL::SubCellspaceInfo* pRPL::SubCellspace::
subInfo() const {
  return _pSubInfo;
}

/*
bool pRPL::SubCellspace::
add2UpdtStream(pRPL::CellStream &updtStream,
               int iDir,
               int iNbr) const {
  if(iDir >= 0 && iDir < _pSubInfo->nNbrDirs() &&
     iNbr >= 0 && iNbr < _pSubInfo->nNbrs(iDir)) {
    if(updtStream.valSize() <= 0) {
      updtStream.init(_pInfo->dataSize());
    }
    else if(updtStream.valSize() != _pInfo->dataSize()) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: Cellspace's data size (" << _pInfo->dataSize() \
           << ") is NOT the same as the CellStream's data size (" \
           << updtStream.valSize() << ")" << endl;
      return false;
    }

    for(int iIdx = 0; iIdx < _vUpdtIdxs.size(); iIdx++) {
      long lclIdx = _vUpdtIdxs[iIdx];
      const pRPL::CoordBR *pSendBR = _pSubInfo->sendBR(iDir, iNbr);
      if(pSendBR != NULL && pSendBR->ifContain(_pInfo->idx2coord(lclIdx))) {
        long glbIdx = _pSubInfo->lclIdx2glbIdx(lclIdx);
        void *aVal = at(lclIdx);
        if(!updtStream.addCell(glbIdx, aVal)) {
          return false;
        }
      }
    } // end -- for(iIdx)
  }
  return true;
}

bool pRPL::SubCellspace::
loadUpdtStream(const pRPL::CellStream &updtStream) {
  if(updtStream.valSize() != _pInfo->dataSize()) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: Cellspace's data size (" << _pInfo->dataSize() \
         << ") is NOT the same as the CellStream's data size (" \
         << updtStream.valSize() << ")" << endl;
    return false;
  }

  long glbIdx, lclIdx;
  for(long iCell = 0; iCell < updtStream.nCells(); iCell++) {
    const void *aCell = updtStream.at(iCell);
    memcpy(&glbIdx, aCell, sizeof(long));
    if(_pSubInfo->MBR().ifContain(_pSubInfo->glbIdx2glbCoord(glbIdx))) {
      lclIdx = _pSubInfo->glbIdx2lclIdx(glbIdx);
      memcpy(at(lclIdx), aCell+sizeof(long), updtStream.valSize());
    }
  }

  return true;
}
*/
