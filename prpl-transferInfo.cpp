#include "prpl-transferInfo.h"

pRPL::TransferInfo::
TransferInfo()
  :_fromPrcID(pRPL::ERROR_ID),
   _toPrcID(pRPL::ERROR_ID),
   _subspcGlbID(pRPL::ERROR_ID),
   _tDataType(pRPL::SPACE_TRANSFDATA),
   _cmplt(false) {}

pRPL::TransferInfo::
TransferInfo(int fromPrcID,
             int toPrcID,
             const string &lyrName,
             int subspcGlbID,
             pRPL::TransferDataType transfType)
  :_fromPrcID(fromPrcID),
   _toPrcID(toPrcID),
   _lyrName(lyrName),
   _subspcGlbID(subspcGlbID),
   _tDataType(transfType),
   _cmplt(false) {}

int pRPL::TransferInfo::
fromPrcID() const {
  return _fromPrcID;
}

int pRPL::TransferInfo::
toPrcID() const {
  return _toPrcID;
}

const string& pRPL::TransferInfo::
lyrName() const {
  return _lyrName;
}

int pRPL::TransferInfo::
subspcGlbID() const {
  return _subspcGlbID;
}

pRPL::TransferDataType pRPL::TransferInfo::
transfDataType() const {
  return _tDataType;
}

bool pRPL::TransferInfo::
completed() const {
  return _cmplt;
}

void pRPL::TransferInfo::
complete() {
  _cmplt = true;
}

void pRPL::TransferInfo::
restart() {
  _cmplt = false;
}

const pRPL::TransferInfo* pRPL::TransferInfoVect::
findInfo(const string &lyrName,
         int subspcGlbID) const {
  const TransferInfo *pInfo = NULL;
  TransferInfoVect::const_iterator itrInfo = begin();
  while(itrInfo != end()) {
    if(itrInfo->lyrName() == lyrName &&
       itrInfo->subspcGlbID() == subspcGlbID) {
      pInfo = &(*itrInfo);
      break;
    }
    itrInfo++;
  }
  return pInfo;
}

bool pRPL::TransferInfoVect::
checkBcastCellspace(int prcID,
                    const string &lyrName) const {
  bool cmplt = true;
  bool found = false;
  for(int iInfo = 0; iInfo < size(); iInfo++) {
    if(at(iInfo).fromPrcID() == prcID &&
       at(iInfo).lyrName() == lyrName &&
       at(iInfo).transfDataType() == pRPL::SPACE_TRANSFDATA &&
       at(iInfo).subspcGlbID() == pRPL::ERROR_ID) {
      found = true;
      if(at(iInfo).completed() == false) {
        cmplt = false;
        break;
      }
    }
  }
  if(found == false) {
    cmplt = false;
  }

  return cmplt;
}
