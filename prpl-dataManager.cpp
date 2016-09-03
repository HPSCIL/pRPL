#include "prpl-dataManager.h"

/****************************************************
*                 Protected Methods                 *
****************************************************/
bool pRPL::DataManager::
_initReadInData(pRPL::Transition &trans,
                pRPL::ReadingOption readOpt) {
  const vector<string> &vInLyrNames = trans.getInLyrNames();
  for(int iInLyr = 0; iInLyr < vInLyrNames.size(); iInLyr++) {
    const string &inLyrName = vInLyrNames.at(iInLyr);
    if(!beginReadingLayer(inLyrName.c_str(), readOpt)) {
      return false;
    }
  } // end -- for(iInLyr) loop

  return true;
}

bool pRPL::DataManager::
_toMsgTag(int &msgTag,
          int lyrID,
          int subspcID,
          int maxSubspcID) {
  if(lyrID < 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__
         << " Error: invalid Layer ID (" << lyrID \
         << ")" << endl;
    return false;
  }
  if((subspcID < 0 && maxSubspcID >= 0) ||
     (subspcID >= 0 && maxSubspcID < 0) ||
     (subspcID >= 0 && maxSubspcID >= 0 && subspcID > maxSubspcID)) {
    cerr << __FILE__ << " function:" << __FUNCTION__
         << " Error: invalid SubCellspace ID (" << subspcID \
         << ") and/or max SubCellspace ID (" << maxSubspcID \
         << ")" << endl;
    return false;
  }

  if(subspcID < 0 && maxSubspcID < 0) { // whole Cellspace
    msgTag = lyrID + 1; // to make sure msgTag is always greater than zero
  }
  else {
    stringstream ssMaxID;
    ssMaxID << maxSubspcID;
    string sMaxID = ssMaxID.str();
    int factor = std::pow(10.0, (int)(sMaxID.size())); // if maxSubspcID = 30, factor = 100;
    msgTag = (lyrID+1) * factor + subspcID;
  }

  return true;
}

bool pRPL::DataManager::
_iBcastCellspaceBegin(const string &lyrName) {
  pRPL::Layer *pLyr = getLayerByName(lyrName.c_str(), true);
  if(pLyr == NULL) {
    _prc.abort();
    return false;
  }

  int msgTag;
  if(!_toMsgTag(msgTag, getLayerIdxByName(lyrName.c_str()),
                pRPL::ERROR_ID, pRPL::ERROR_ID)) {
    _prc.abort();
    return false;
  }

  if(_prc.isMaster()) { // Master process
    pRPL::Cellspace *pSpc = pLyr->cellspace();
    if(pSpc == NULL || pSpc->isEmpty()) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: NULL or empty Cellspace on Process ("
          << _prc.id() << ") to transfer" \
          << endl;
      _prc.abort();
      return false;
    }

    int msgSize = pSpc->info()->size() * pSpc->info()->dataSize();
    for(int toPrcID = 1; toPrcID < _prc.nProcesses(); toPrcID++) {
      _vBcstSpcReqs.resize(_vBcstSpcReqs.size() + 1);
      MPI_Isend(pSpc->_pData, msgSize, MPI_CHAR, toPrcID, msgTag, _prc.comm(), &(_vBcstSpcReqs.back()));
      _vBcstSpcInfos.push_back(pRPL::TransferInfo(_prc.id(), toPrcID,
                                                  lyrName, pRPL::ERROR_ID,
                                                  pRPL::SPACE_TRANSFDATA));
    }
  } // end -- if(_prc.isMaster())

  else if(!_prc.isWriter()) { // Worker process
    pRPL::Cellspace *pSpc = pLyr->cellspace();
    if(pSpc == NULL) {
      pLyr->createCellspace();
      pSpc = pLyr->cellspace();
    }
    if(pSpc == NULL || pSpc->isEmpty()) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: NULL or empty Cellspace on Process ("
          << _prc.id() << ") to transfer" \
          << endl;
      _prc.abort();
      return false;
    }

    int msgSize = pSpc->info()->size() * pSpc->info()->dataSize();
    _vBcstSpcReqs.resize(_vBcstSpcReqs.size() + 1);
    MPI_Irecv(pSpc->_pData, msgSize, MPI_CHAR, 0, msgTag, _prc.comm(), &(_vRecvSpcReqs.back()));
    _vBcstSpcInfos.push_back(pRPL::TransferInfo(0, _prc.id(),
                                                lyrName, pRPL::ERROR_ID,
                                                pRPL::SPACE_TRANSFDATA));
  } // end -- worker process

  return true;
}

bool pRPL::DataManager::
_iTransfSubspaceBegin(int fromPrcID,
                      int toPrcID,
                      const string &lyrName,
                      int subspcGlbID) {
  pRPL::Layer *pLyr = getLayerByName(lyrName.c_str(), true);
  if(pLyr == NULL) {
    _prc.abort();
    return false;
  }

  int maxSubspcID = pLyr->maxSubCellspaceID();
  if(subspcGlbID < 0 || subspcGlbID > maxSubspcID) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: invalid SubCellspace ID (" << subspcGlbID \
        << ")" << endl;
    _prc.abort();
    return false;
  }

  int msgTag;
  if(!_toMsgTag(msgTag, getLayerIdxByName(lyrName.c_str()),
                subspcGlbID, maxSubspcID)) {
    _prc.abort();
    return false;
  }

  if(_prc.id() == fromPrcID) { // FROM process
    pRPL::Cellspace *pSpc = pLyr->subCellspace_glbID(subspcGlbID, true);
    if(pSpc == NULL || pSpc->isEmpty()) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: NULL or empty SubCellspace on Process ("
           << _prc.id() << ") to transfer" \
           << endl;
      _prc.abort();
      return false;
    }

    int msgSize = pSpc->info()->size() * pSpc->info()->dataSize();
    _vSendSpcReqs.resize(_vSendSpcReqs.size() + 1);
    MPI_Isend(pSpc->_pData, msgSize, MPI_CHAR, toPrcID, msgTag, _prc.comm(), &(_vSendSpcReqs.back()));
    _vSendSpcInfos.push_back(pRPL::TransferInfo(_prc.id(), toPrcID,
                             lyrName, subspcGlbID,
                             pRPL::SPACE_TRANSFDATA));
  } // end -- if(_prc.id() == fromPrcID)

  else if(_prc.id() == toPrcID) { // TO process
    pRPL::Cellspace *pSpc = pLyr->subCellspace_glbID(subspcGlbID);
    if(pSpc == NULL) {
      pSpc = pLyr->addSubCellspace(subspcGlbID);
    }
    if(pSpc == NULL || pSpc->isEmpty()) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: NULL or empty SubCellspace on Process ("
          << _prc.id() << ") to transfer" \
          << endl;
      _prc.abort();
      return false;
    }

    int msgSize = pSpc->info()->size() * pSpc->info()->dataSize();
    _vRecvSpcReqs.resize(_vRecvSpcReqs.size() + 1);
    MPI_Irecv(pSpc->_pData, msgSize, MPI_CHAR, fromPrcID, msgTag, _prc.comm(), &(_vRecvSpcReqs.back()));
    _vRecvSpcInfos.push_back(pRPL::TransferInfo(fromPrcID, _prc.id(),
                                                lyrName, subspcGlbID,
                                                pRPL::SPACE_TRANSFDATA));
  } // end -- if(_prc.id() == toPrcID)

  return true;
}

int pRPL::DataManager::
_iTransfSpaceTest(pRPL::IntVect &vCmpltIDs,
                  pRPL::TransferType transfType,
                  bool delSentSpaces,
                  bool writeRecvSpaces,
                  bool delRecvSpaces) {
  int nCmplts = 0;

  vector<MPI_Request> *pvTransfReqs = NULL;
  pRPL::TransferInfoVect *pvTransfInfos = NULL;
  switch(transfType) {
    case pRPL::BCST_TRANSF:
      pvTransfReqs = &_vBcstSpcReqs;
      pvTransfInfos = &_vBcstSpcInfos;
      break;
    case pRPL::SEND_TRANSF:
      pvTransfReqs = &_vSendSpcReqs;
      pvTransfInfos = &_vSendSpcInfos;
      break;
    case pRPL::RECV_TRANSF:
      pvTransfReqs = &_vRecvSpcReqs;
      pvTransfInfos = &_vRecvSpcInfos;
      break;
    default:
      cerr << __FILE__ << " function:" << __FUNCTION__
           << " Error: invalid TransferType (" << transfType \
           << ")" << endl;
      return MPI_UNDEFINED;
  }

  if(pvTransfReqs->empty()) {
    return MPI_UNDEFINED;
  }

  vCmpltIDs.clear();
  vCmpltIDs.resize(pvTransfReqs->size());
  vector<MPI_Status> vCmpltStats(pvTransfReqs->size());

  MPI_Testsome(pvTransfReqs->size(), &(pvTransfReqs->at(0)),
               &nCmplts, &(vCmpltIDs[0]), &(vCmpltStats[0]));
  if(nCmplts >= 0) {
    vCmpltIDs.resize(nCmplts);
    vCmpltStats.resize(nCmplts);
  }

  for(int iCmplt = 0; iCmplt < nCmplts; iCmplt++) {
    pRPL::TransferInfo &transfInfo = pvTransfInfos->at(vCmpltIDs[iCmplt]);
    transfInfo.complete();

    int fromPrcID = transfInfo.fromPrcID();
    int toPrcID = transfInfo.toPrcID();
    int subspcGlbID = transfInfo.subspcGlbID();

    const string &lyrName = transfInfo.lyrName();
    pRPL::Layer *pLyr = getLayerByName(lyrName.c_str(), true);
    if(pLyr == NULL) {
      continue;
    }

    if(delSentSpaces &&
       _prc.id() == fromPrcID) {
      if(transfType == pRPL::SEND_TRANSF) {
        if(!pLyr->delSubCellspace_glbID(subspcGlbID)){
          cerr << __FILE__ << " function:" << __FUNCTION__
               << " Error: unable to find SubCellspace (" \
               << subspcGlbID << ") to delete" << endl;
        }
      }
      else if(transfType == pRPL::BCST_TRANSF) {
        if(pvTransfInfos->checkBcastCellspace(_prc.id(), lyrName)) {
          pLyr->delCellspace();
        }
      }
    } // end -- if(delSentSpaces)

    if(writeRecvSpaces &&
       _prc.id() == toPrcID &&
       transfType == pRPL::RECV_TRANSF) {
      //cout << _prc.id() << ": writing SubCellspace " << subspcGlbID << endl;
      if(!pLyr->writeSubCellspaceByGDAL(subspcGlbID)) {
        cerr << __FILE__ << " function:" << __FUNCTION__
             << " Error: failed to write SubCellspace (" \
             << subspcGlbID << ")" << endl;
      }
    } // end -- if(writeRecvSpaces)

    if(delRecvSpaces &&
       _prc.id() == toPrcID &&
       transfType == pRPL::RECV_TRANSF) {
      if(!pLyr->delSubCellspace_glbID(subspcGlbID)) {
        cerr << __FILE__ << " function:" << __FUNCTION__
             << " Error: unable to find SubCellspace (" \
             << subspcGlbID << ") to delete" << endl;
      }
    } // end -- if(delRecvSpaces)

  } // end -- for(iCmplt) loop

  return nCmplts;
}

void pRPL::DataManager::
_clearTransfSpaceReqs() {
  _vBcstSpcReqs.clear();
  _vBcstSpcInfos.clear();
  _vSendSpcReqs.clear();
  _vSendSpcInfos.clear();
  _vRecvSpcReqs.clear();
  _vRecvSpcInfos.clear();
  _vEvaledSubspcIDs.clear();
}

int pRPL::DataManager::
_checkTransInDataReady(pRPL::Transition &trans,
                       pRPL::ReadingOption readOpt) {
  int readySubspcID = pRPL::ERROR_ID, subspcID2Check;

  const pRPL::IntVect *pvSubspcIDs = _mOwnerships.subspcIDsOnPrc(_prc.id());
  if(pvSubspcIDs == NULL) {
    return pRPL::ERROR_ID;
  }

  for(int iSubspc = 0; iSubspc < pvSubspcIDs->size(); iSubspc++) {
    subspcID2Check = pvSubspcIDs->at(iSubspc);
    if(std::find(_vEvaledSubspcIDs.begin(), _vEvaledSubspcIDs.end(), subspcID2Check) != _vEvaledSubspcIDs.end()) {
      continue; // this SubCellspace has been evaluated, ignore the following steps, check the next SubCellspace
    }
    else { // if SubCellspace has NOT been evaluated
      if(readOpt == pRPL::PARA_READING) { // if parallel reading, input data has been read
        readySubspcID = subspcID2Check;
        break; // break the for(iSubspc) loop
      }
      else if(readOpt == pRPL::CENT_READING ||
              readOpt == pRPL::CENTDEL_READING) { // if centralized reading, check if transfers have completed
        const vector<string> &vInLyrNames = trans.getInLyrNames();
        for(int iInLyr = 0; iInLyr < vInLyrNames.size(); iInLyr++) {
          pRPL::Layer *pInLyr = getLayerByName(vInLyrNames[iInLyr].c_str(), true);
          if(pInLyr == NULL) {
            return pRPL::ERROR_ID;
          }
          if(pInLyr->isDecomposed()) {
            const pRPL::TransferInfo *pTransfInfo = _vRecvSpcInfos.findInfo(vInLyrNames[iInLyr], subspcID2Check);
            if(pTransfInfo == NULL || pTransfInfo->completed() == false) {
              subspcID2Check = pRPL::ERROR_ID;
              break; // the transfer of this SubCellspace has NOT been completed, break the for(iInLyr) loop
            }
          }
          else {
            const pRPL::TransferInfo *pTransfInfo = _vBcstSpcInfos.findInfo(vInLyrNames[iInLyr], pRPL::ERROR_ID);
            if(pTransfInfo == NULL || pTransfInfo->completed() == false) {
              subspcID2Check = pRPL::ERROR_ID;
              break; // the transfer (broadcast) of this Cellspace has NOT been completed, break the for(iInLyr) loop
            }
          }
        } // end -- for(iInLyr) loop

        if(subspcID2Check != pRPL::ERROR_ID) {
          readySubspcID = subspcID2Check;
          break; // found a SubCellspace whose transfer has been completed, break the for(iSubspc) loop
        }
      } // end -- if(readOpt == pRPL::CENT_READING)
      else if(readOpt == pRPL::PGT_READING) {
        readySubspcID = subspcID2Check;
        break; // break the for(iSubspc) loop
      }
    } // end -- if SubCellspace has NOT been evaluated
  } // end -- for(iSubspc) loop

  return readySubspcID;
}

bool pRPL::DataManager::
_initTransOutData(pRPL::Transition &trans,
                  int subspcGlbID,
                  bool ifInitNoData) {
  const vector<string> &vOutLyrNames = trans.getOutLyrNames();
  for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
    pRPL::Layer *pOutLyr = getLayerByName(vOutLyrNames[iOutLyr].c_str(), true);
    if(pOutLyr == NULL) {
      return false;
    }
    if(pOutLyr->cellspaceInfo() == NULL) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: Layer (" << vOutLyrNames[iOutLyr] \
           << ") does NOT have a global CellspaceInfo initialized" \
           << endl;
      return false;
    }

    double *pNoDataVal = NULL;
    double noDataVal;
    if(ifInitNoData) {
      noDataVal = pOutLyr->cellspaceInfo()->getNoDataValAs<double>();
      pNoDataVal = &noDataVal;
    }

    if(pOutLyr->isDecomposed()) { // decomposed Layer
      pRPL::SubCellspace *pSubspc = pOutLyr->subCellspace_glbID(subspcGlbID, false);
      if(pSubspc == NULL) {
        pSubspc = pOutLyr->addSubCellspace(subspcGlbID, pNoDataVal);
      }
      if(pSubspc == NULL) {
        return false;
      }
    }
    else { // whole Cellspace on the Layer
      pRPL::Cellspace *pSpc = pOutLyr->cellspace();
      if(pSpc == NULL) {
        pOutLyr->createCellspace(pNoDataVal);
        pSpc = pOutLyr->cellspace();
      }
      if(pSpc == NULL) {
        return false;
      }
    }
  } // end -- for(iOutLyr)

  return true;
}

bool pRPL::DataManager::
_setTransCellspaces(pRPL::Transition &trans,
                    int subspcGlbID) {
  trans.clearCellspaces();

  const vector<string> &vInLyrNames = trans.getInLyrNames();
  for(size_t iLyr = 0; iLyr < vInLyrNames.size(); iLyr++) {
    pRPL::Layer *pLyr = getLayerByName(vInLyrNames[iLyr].c_str(), true);
    if(pLyr == NULL) {
      return false;
    }
    if(!pLyr->isDecomposed() &&
       pLyr->cellspace() != NULL &&
       !pLyr->cellspace()->isEmpty()) { // use whole Cellspace
      if(!trans.setCellspace(vInLyrNames[iLyr], pLyr->cellspace())) {
        return false;
      }
    }
    else { // use SubCellspace
      pRPL::SubCellspace *pSubSpc = pLyr->subCellspace_glbID(subspcGlbID, true);
      if(pSubSpc == NULL) {
        return false;
      }
      if(!trans.setCellspace(vInLyrNames[iLyr], pSubSpc)) {
        return false;
      }
    }
  } // end -- for(iLyr)

  const vector<string> &vOutLyrNames = trans.getOutLyrNames();
  for(size_t iLyr = 0; iLyr < vOutLyrNames.size(); iLyr++) {
    if(trans.getCellspaceByLyrName(vOutLyrNames[iLyr]) != NULL) {
      continue; // ignore the SubCellspaces that have been set as input data
    }
    pRPL::Layer *pLyr = getLayerByName(vOutLyrNames[iLyr].c_str(), true);
    if(pLyr == NULL) {
      return false;
    }
    if(!pLyr->isDecomposed() &&
       pLyr->cellspace() != NULL &&
       !pLyr->cellspace()->isEmpty()) { // use whole Cellspace
      if(!trans.setCellspace(vOutLyrNames[iLyr], pLyr->cellspace())) {
        return false;
      }
    }
    else { // use SubCellspace
      SubCellspace *pSubSpc = pLyr->subCellspace_glbID(subspcGlbID, true);
      if(pSubSpc == NULL) {
        return false;
      }
      if(!trans.setCellspace(vOutLyrNames[iLyr], pSubSpc)) {
        return false;
      }
    }
  } // end -- for(iLyr)

  if(!trans.afterSetCellspaces(subspcGlbID)) {
    return false;
  }

  return true;
}

bool pRPL::DataManager::
_setTransNbrhd(pRPL::Transition &trans) {
  trans.setNbrhd();

  if(!trans.getNbrhdName().empty()) {
    pRPL::Neighborhood *pNbrhd = getNbrhdByName(trans.getNbrhdName().c_str(), true);
    if(pNbrhd != NULL) {
      trans.setNbrhd(pNbrhd);
      if(!trans.afterSetNbrhd()) {
        return false;
      }
    }
    else {
      return false;
    }
  }

  return true;
}

bool pRPL::DataManager::
_initEvaluate(pRPL::Transition &trans,
              int subspcGlbID,
              pRPL::EvaluateBR br2Eval,
              bool ifInitNoData) {
  if(!_setTransNbrhd(trans)) {
    return false;
  }
  if(br2Eval != pRPL::EVAL_INTERIOR) {
    if(!_initTransOutData(trans, subspcGlbID, ifInitNoData)) {
      return false;
    }
  }
  if(!_setTransCellspaces(trans, subspcGlbID)) {
    return false;
  }

  if(!trans.check()) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: Transition check failed" << endl;
    return false;
  }

  if(br2Eval != pRPL::EVAL_INTERIOR && trans.needExchange()) {
    trans.setUpdateTracking(true);
  }

  return true;
}

bool pRPL::DataManager::
_fnlzEvaluate(pRPL::Transition &trans,
              int subspcGlbID) {
  trans.setUpdateTracking(false);
  return true;
}

bool pRPL::DataManager::
_workerReadSubspcs(pRPL::Transition &trans,
                   int subspcGlbID,
                   pRPL::ReadingOption readOpt) {
  const vector<string> &vInLyrNames = trans.getInLyrNames();
  for(int iLyr = 0; iLyr < vInLyrNames.size(); iLyr++) {
    const string &inLyrName = vInLyrNames[iLyr];
    pRPL::Layer *pLyr = getLayerByName(inLyrName.c_str(), true);
    if(pLyr == NULL) {
      return false;
    }

    if(pLyr->isDecomposed()) { // decomposed Layer
      if(readOpt == pRPL::PARA_READING) { // parallel reading
        if(!pLyr->loadSubCellspaceByGDAL(subspcGlbID)) {
          return false;
        }
      }
      else if(readOpt == pRPL::CENT_READING ||
              readOpt == pRPL::CENTDEL_READING) { // centralized reading
        if(!_iTransfSubspaceBegin(0, _prc.id(), inLyrName, subspcGlbID)) {
          return false;
        }
      }
      else if(readOpt == pRPL::PGT_READING) {
        if(!pLyr->loadSubCellspaceByPGTIOL(subspcGlbID)) {
          return false;
        }
      }
    }
    else { // non-decomposed Layer
      cerr << __FILE__ << " function:" << __FUNCTION__
           << " Error: unable to read a SubCellspace within a non-decomposed Layer (" \
           << inLyrName << ")" \
           << endl;
      return false;
    }
  } // end -- for(iLyr)
  return true;
}

bool pRPL::DataManager::
_workerWriteSubspcs(pRPL::Transition &trans,
                    int subspcGlbID,
                    pRPL::WritingOption writeOpt) {
  if(writeOpt == pRPL::NO_WRITING) {
    return true;
  }

  pRPL::TaskSignal req = pRPL::SUBMIT_SGNL;
  int submit2ID = _prc.hasWriter() ? 1 : 0; // if a writer exists, send to writer; otherwise master
  const vector<string> &vOutLyrNames = trans.getOutLyrNames();
  for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
    const string &outLyrName = vOutLyrNames[iOutLyr];
    pRPL::Layer *pOutLyr = getLayerByName(outLyrName.c_str(), true);
    if(pOutLyr == NULL) {
      return false;
    }
    if(!pOutLyr->isDecomposed()) {
      cerr << __FILE__ << " function:" << __FUNCTION__
           << " Error: unable to write a SubCellspace within a non-decomposed Layer (" \
           << outLyrName << ")" \
           << endl;
      return false;
    }

    if(writeOpt == pRPL::PARA_WRITING ||
       writeOpt == pRPL::PARADEL_WRITING) {
      if(!pOutLyr->writeTmpSubCellspaceByGDAL(subspcGlbID, _prc.nProcesses())) {
        return false;
      }
      if(writeOpt == pRPL::PARADEL_WRITING) {
        if(!pOutLyr->delSubCellspace_glbID(subspcGlbID)){
          cerr << __FILE__ << " function:" << __FUNCTION__
               << " Error: unable to find SubCellspace (" \
               << subspcGlbID << ") to delete" << endl;
        }
      }
    }
    else if(writeOpt == pRPL::CENT_WRITING ||
            writeOpt == pRPL::CENTDEL_WRITING) {
      int lyrID = getLayerIdxByName(outLyrName.c_str(), true);
      MPI_Send(&req, 1, MPI_INT, submit2ID, pRPL::REQUEST_TAG, _prc.comm());
      int sendSubspcInfo[2];
      sendSubspcInfo[0] = lyrID;
      sendSubspcInfo[1] = subspcGlbID;
      MPI_Send(sendSubspcInfo, 2, MPI_INT, submit2ID, pRPL::SUBMIT_TAG, _prc.comm());
      if(!_iTransfSubspaceBegin(_prc.id(), submit2ID, outLyrName, subspcGlbID)) {
        return false;
      }
    }
    else if(writeOpt == pRPL::PGT_WRITING ||
            writeOpt == pRPL::PGTDEL_WRITING) {
      if(!pOutLyr->writeSubCellspaceByPGTIOL(subspcGlbID)) {
        return false;
      }

      if(writeOpt == pRPL::PGTDEL_WRITING) {
        if(!pOutLyr->delSubCellspace_glbID(subspcGlbID)){
          cerr << __FILE__ << " function:" << __FUNCTION__
               << " Error: unable to find SubCellspace (" \
               << subspcGlbID << ") to delete" << endl;
        }
      }
    }
  } // end -- for(iOutLyr)

  if(writeOpt == pRPL::PARA_WRITING ||
     writeOpt == pRPL::PARADEL_WRITING ||
     writeOpt == pRPL::PGT_WRITING ||
     writeOpt == pRPL::PGTDEL_WRITING) {
    MPI_Send(&req, 1, MPI_INT, submit2ID, pRPL::REQUEST_TAG, _prc.comm());
    MPI_Send(&subspcGlbID, 1, MPI_INT, submit2ID, pRPL::SUBMIT_TAG, _prc.comm());
  }

  return true;
}

bool pRPL::DataManager::
_dymWriteSubspcs(pRPL::Transition &trans,
                 int subspcGlbID,
                 pRPL::WritingOption writeOpt) {
  if(writeOpt == pRPL::NO_WRITING) {
    return true;
  }

  pRPL::TaskSignal req = pRPL::SUBMIT_SGNL;
  const vector<string> &vOutLyrNames = trans.getOutLyrNames();
  for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
    const string &outLyrName = vOutLyrNames[iOutLyr];
    pRPL::Layer *pOutLyr = getLayerByName(outLyrName.c_str(), true);
    if(pOutLyr == NULL) {
      return false;
    }
    if(!pOutLyr->isDecomposed()) {
      cerr << __FILE__ << " function:" << __FUNCTION__
          << " Error: unable to write a SubCellspace within a non-decomposed Layer (" \
          << outLyrName << ")" \
          << endl;
      return false;
    }

    if(writeOpt == pRPL::PARA_WRITING ||
       writeOpt == pRPL::PARADEL_WRITING) {
      if(!pOutLyr->writeTmpSubCellspaceByGDAL(subspcGlbID, _prc.nProcesses())) {
        return false;
      }
      if(writeOpt == pRPL::PARADEL_WRITING) {
        if(!pOutLyr->delSubCellspace_glbID(subspcGlbID)){
          cerr << __FILE__ << " function:" << __FUNCTION__
               << " Error: unable to find SubCellspace (" \
               << subspcGlbID << ") to delete" << endl;
        }
      }
    }
    else if((writeOpt == pRPL::CENT_WRITING || writeOpt == pRPL::CENTDEL_WRITING) &&
            _prc.hasWriter()) {
      int lyrID = getLayerIdxByName(outLyrName.c_str(), true);
      MPI_Send(&req, 1, MPI_INT, 1, pRPL::REQUEST_TAG, _prc.comm());
      int sendSubspcInfo[2];
      sendSubspcInfo[0] = lyrID;
      sendSubspcInfo[1] = subspcGlbID;
      MPI_Send(sendSubspcInfo, 2, MPI_INT, 1, pRPL::SUBMIT_TAG, _prc.comm());
      if(!_iTransfSubspaceBegin(_prc.id(), 1, outLyrName, subspcGlbID)) {
        return false;
      }
    }
    else if(writeOpt == pRPL::PGT_WRITING ||
            writeOpt == pRPL::PGTDEL_WRITING) {
      if(!pOutLyr->writeSubCellspaceByPGTIOL(subspcGlbID)) {
        return false;
      }
      if(writeOpt == pRPL::PGTDEL_WRITING) {
        if(!pOutLyr->delSubCellspace_glbID(subspcGlbID)) {
          cerr << __FILE__ << " function:" << __FUNCTION__
              << " Error: unable to find SubCellspace (" \
              << subspcGlbID << ") to delete" << endl;
        }
      }
    }
  } // end -- for(iOutLyr)

  if((writeOpt == pRPL::PARA_WRITING ||
      writeOpt == pRPL::PARADEL_WRITING ||
      writeOpt == pRPL::PGT_WRITING ||
      writeOpt == pRPL::PGTDEL_WRITING) &&
     _prc.hasWriter()) {
    MPI_Send(&req, 1, MPI_INT, 1, pRPL::REQUEST_TAG, _prc.comm());
    MPI_Send(&subspcGlbID, 1, MPI_INT, 1, pRPL::SUBMIT_TAG, _prc.comm());
  }

  return true;
}

bool pRPL::DataManager::
_finalWriteSubspcs(pRPL::Transition &trans,
                   pRPL::WritingOption writeOpt) {
  if(writeOpt == pRPL::NO_WRITING) {
    return true;
  }

  if((writeOpt == pRPL::CENT_WRITING || writeOpt == pRPL::CENTDEL_WRITING) &&
     !_prc.hasWriter() &&
     !_prc.isMaster()) {
    pRPL::TaskSignal req = pRPL::SUBMIT_SGNL;
    const IntVect *pvSubspcIDs = _mOwnerships.subspcIDsOnPrc(_prc.id());
    if(pvSubspcIDs != NULL) {
      for(int iSubspc = 0; iSubspc < pvSubspcIDs->size(); iSubspc++) {
        int subspcGlbID = pvSubspcIDs->at(iSubspc);

        const vector<string> &vOutLyrNames = trans.getOutLyrNames();
        for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
          const string &outLyrName = vOutLyrNames[iOutLyr];
          pRPL::Layer *pOutLyr = getLayerByName(outLyrName.c_str(), true);
          if(pOutLyr == NULL) {
            return false;
          }
          if(!pOutLyr->isDecomposed()) {
            cerr << __FILE__ << " function:" << __FUNCTION__
                << " Error: unable to write a SubCellspace within a non-decomposed Layer (" \
                << outLyrName << ")" \
                << endl;
            return false;
          }

          int lyrID = getLayerIdxByName(outLyrName.c_str(), true);
          MPI_Send(&req, 1, MPI_INT, 0, pRPL::REQUEST_TAG, _prc.comm());
          int sendSubspcInfo[2];
          sendSubspcInfo[0] = lyrID;
          sendSubspcInfo[1] = subspcGlbID;
          MPI_Send(sendSubspcInfo, 2, MPI_INT, 0, pRPL::SUBMIT_TAG, _prc.comm());
          if(!_iTransfSubspaceBegin(_prc.id(), 0, outLyrName, subspcGlbID)) {
            return false;
          }

        } // end -- for(iOutLyr)
      } // end -- for(iSubspc)
    } // end -- if(pvSubspcIDs != NULL)
  }

  return true;
}

pRPL::EvaluateReturn pRPL::DataManager::
_evalNone(pRPL::Transition &trans,
          int subspcGlbID,
          bool ifInitNoData) {
  if(!_initEvaluate(trans, subspcGlbID, pRPL::EVAL_WORKBR, ifInitNoData)) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to initialize the Transition" \
         << endl;
    return pRPL::EVAL_FAILED;
  }
  if(!_fnlzEvaluate(trans, subspcGlbID)) {
    return pRPL::EVAL_FAILED;
  }

  return pRPL::EVAL_SUCCEEDED;
}

pRPL::EvaluateReturn pRPL::DataManager::
_evalAll(pRPL::Transition &trans,
         int subspcGlbID,
         pRPL::EvaluateBR br2Eval,
         bool ifInitNoData) {
  if(br2Eval != pRPL::EVAL_WORKBR) {
    if(subspcGlbID == pRPL::ERROR_ID) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: cannot evaluate the Edges or Interior of a whole Cellspace" \
           << endl;
      return pRPL::EVAL_FAILED;
    }
    if(trans.getOutLyrNames().empty()) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: no output Layers to evaluate Edges or Interior" \
           << endl;
      return pRPL::EVAL_FAILED;
    }
    if(!trans.isOutLyr(trans.getPrimeLyrName())) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: the primary Layer (" << trans.getPrimeLyrName() \
           << ") is NOT an output Layer" << endl;
      return pRPL::EVAL_FAILED;
    }
  }

  if(!_initEvaluate(trans, subspcGlbID, br2Eval, ifInitNoData)) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to initialize the Transition" \
        << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::Layer *pPrmLyr = getLayerByName(trans.getPrimeLyrName().c_str(), true);
  if(pPrmLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to find the primary Layer" \
         << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::Cellspace *pPrmCellspace = NULL;
  if(subspcGlbID == pRPL::ERROR_ID) { // if evaluate the whole Cellspace
    pPrmCellspace = pPrmLyr->cellspace();
  }
  else { // if evaluate a SubCellspace
    pPrmCellspace = pPrmLyr->subCellspace_glbID(subspcGlbID, true);
  }
  if(pPrmCellspace == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to find the primary Cellspace" \
         << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::CoordBR workBR;
  const pRPL::CoordBR *pWorkBR;
  pRPL::EvaluateReturn done = pRPL::EVAL_SUCCEEDED;
  switch(br2Eval) {
    case pRPL::EVAL_WORKBR:
      if(subspcGlbID == pRPL::ERROR_ID) { // if evaluate the whole Cellspace
        if(!pPrmCellspace->info()->calcWorkBR(&workBR, trans.getNbrhd(), true)) {
          return pRPL::EVAL_FAILED;
        }
      }
      else { // if evaluate a SubCellspace
        workBR = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->workBR();
      }
      done = trans.evalBR(workBR);
      break;
    case pRPL::EVAL_EDGES:
      for(int iEdge = 0; iEdge < ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->nEdges(); iEdge++) {
        pWorkBR = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->edgeBR(iEdge);
        if(pWorkBR != NULL) {
          done = trans.evalBR(*pWorkBR);
          if(done == pRPL::EVAL_FAILED) {
            break;
          }
        }
      }
      break;
    case pRPL::EVAL_INTERIOR:
      pWorkBR = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->interiorBR();
      if(pWorkBR != NULL) {
        done = trans.evalBR(*pWorkBR);
      }
      break;
    default:
      break;
  }

  if(!_fnlzEvaluate(trans, subspcGlbID)) {
    return pRPL::EVAL_FAILED;
  }

  return done;
}

pRPL::EvaluateReturn pRPL::DataManager::
_evalRandomly(pRPL::Transition &trans,
              int subspcGlbID,
              bool ifInitNoData) {
  if(!_initEvaluate(trans, subspcGlbID, pRPL::EVAL_WORKBR, ifInitNoData)) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to initialize the Transition" \
         << endl;
    return pRPL::EVAL_FAILED;
  }
  if(trans.needExchange()) {
    trans.setUpdateTracking(true);
  }

  pRPL::Layer *pPrmLyr = getLayerByName(trans.getPrimeLyrName().c_str(), true);
  if(pPrmLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to find the primary Layer" \
         << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::Cellspace *pPrmCellspace = NULL;
  if(subspcGlbID == pRPL::ERROR_ID) { // if evaluate the whole Cellspace
    pPrmCellspace = pPrmLyr->cellspace();
  }
  else { // if evaluate a SubCellspace
    pPrmCellspace = pPrmLyr->subCellspace_glbID(subspcGlbID, true);
  }
  if(pPrmCellspace == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to find the primary Cellspace" \
        << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::CoordBR workBR;
  if(subspcGlbID == pRPL::ERROR_ID) { // if evaluate the whole Cellspace
    if(!pPrmCellspace->info()->calcWorkBR(&workBR, trans.getNbrhd(), true)) {
      return pRPL::EVAL_FAILED;
    }
  }
  else { // if evaluate a SubCellspace
    workBR = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->workBR();
  }
  pRPL::EvaluateReturn done = trans.evalRandomly(workBR);

  if(!_fnlzEvaluate(trans, subspcGlbID)) {
    return pRPL::EVAL_FAILED;
  }

  return done;
}

pRPL::EvaluateReturn pRPL::DataManager::
_evalSelected(pRPL::Transition &trans,
              const pRPL::LongVect &vGlbIdxs,
              int subspcGlbID,
              pRPL::EvaluateBR br2Eval,
              bool ifInitNoData) {
  if(br2Eval != pRPL::EVAL_WORKBR) {
    if(subspcGlbID == pRPL::ERROR_ID) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: cannot evaluate the Edges or Interior of a whole Cellspace" \
          << endl;
      return pRPL::EVAL_FAILED;
    }
    if(trans.getOutLyrNames().empty()) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: no output Layers to evaluate Edges or Interior" \
           << endl;
      return pRPL::EVAL_FAILED;
    }
    if(!trans.isOutLyr(trans.getPrimeLyrName())) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: the primary Layer (" << trans.getPrimeLyrName() \
           << ") is NOT an output Layer" << endl;
      return pRPL::EVAL_FAILED;
    }
  }

  if(!_initEvaluate(trans, subspcGlbID, br2Eval, ifInitNoData)) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to initialize the Transition" \
        << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::Layer *pPrmLyr = getLayerByName(trans.getPrimeLyrName().c_str(), true);
  if(pPrmLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to find the primary Layer" \
        << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::Cellspace *pPrmCellspace = NULL;
  if(subspcGlbID == pRPL::ERROR_ID) { // if evaluate the whole Cellspace
    pPrmCellspace = pPrmLyr->cellspace();
  }
  else { // if evaluate a SubCellspace
    pPrmCellspace = pPrmLyr->subCellspace_glbID(subspcGlbID, true);
  }
  if(pPrmCellspace == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to find the primary Cellspace" \
        << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::CoordBR workBR;
  const pRPL::CoordBR *pWorkBR;
  pRPL::EvaluateReturn done = pRPL::EVAL_SUCCEEDED;
  switch(br2Eval) {
    case pRPL::EVAL_WORKBR:
      if(subspcGlbID == pRPL::ERROR_ID) { // if evaluate the whole Cellspace
        if(!pPrmCellspace->info()->calcWorkBR(&workBR, trans.getNbrhd(), true)) {
          return pRPL::EVAL_FAILED;
        }
        done = trans.evalSelected(workBR, vGlbIdxs);
      }
      else { // if evaluate a SubCellspace
        pWorkBR = &(((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->workBR());
        for(int iIdx = 0; iIdx < vGlbIdxs.size(); iIdx++) {
          pRPL::CellCoord lclCoord = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->glbIdx2lclCoord(vGlbIdxs[iIdx]);
          if(pWorkBR->ifContain(lclCoord)) {
            done = trans.evaluate(lclCoord);
            if(done == pRPL::EVAL_FAILED ||
               done == pRPL::EVAL_TERMINATED) {
              return done;
            }
          }
        }
      }
      break;
    case pRPL::EVAL_EDGES:
      for(int iEdge = 0; iEdge < ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->nEdges(); iEdge++) {
        pWorkBR = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->edgeBR(iEdge);
        if(pWorkBR != NULL) {
          for(int iIdx = 0; iIdx < vGlbIdxs.size(); iIdx++) {
            pRPL::CellCoord lclCoord = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->glbIdx2lclCoord(vGlbIdxs[iIdx]);
            if(pWorkBR->ifContain(lclCoord)) {
              done = trans.evaluate(lclCoord);
              if(done == pRPL::EVAL_FAILED ||
                 done == pRPL::EVAL_TERMINATED) {
                return done;
              }
            }
          } // end -- for(iIdx)
        }
      } // end -- for(iEdge)
      break;
    case pRPL::EVAL_INTERIOR:
      pWorkBR = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->interiorBR();
      if(pWorkBR != NULL) {
        for(int iIdx = 0; iIdx < vGlbIdxs.size(); iIdx++) {
          pRPL::CellCoord lclCoord = ((pRPL::SubCellspace *)pPrmCellspace)->subInfo()->glbIdx2lclCoord(vGlbIdxs[iIdx]);
          if(pWorkBR->ifContain(lclCoord)) {
            done = trans.evaluate(lclCoord);
            if(done == pRPL::EVAL_FAILED ||
               done == pRPL::EVAL_TERMINATED) {
              return done;
            }
          }
        }
      }
      break;
    default:
      break;
  }

  if(!_fnlzEvaluate(trans, subspcGlbID)) {
    return pRPL::EVAL_FAILED;
  }

  return done;
}

pRPL::EvaluateReturn pRPL::DataManager::
_masterTF(pRPL::Transition &trans,
		  pRPL::ReadingOption readOpt,
		  pRPL::WritingOption writeOpt) {
  int nCmpltBcsts = 0, nCmpltSends = 0, nCmpltRecvs = 0;
  pRPL::IntVect vCmpltBcstIDs, vCmpltSendIDs, vCmpltRecvIDs;
  pRPL::TaskSignal req;
  MPI_Status status;

  bool delSentSpc = (readOpt == pRPL::CENTDEL_READING) ? true : false;
  bool writeRecvSpc = (writeOpt != pRPL::NO_WRITING) ? true : false;
  bool delRecvSpc = (writeOpt == pRPL::CENTDEL_WRITING) ? true : false;

  while(_nActiveWorkers > 0) {
    MPI_Recv(&req, 1, MPI_INT, MPI_ANY_SOURCE, pRPL::REQUEST_TAG, _prc.comm(), &status);
    int workerID = status.MPI_SOURCE;

    if((pRPL::TaskSignal)req == pRPL::NEWTASK_SGNL) { // worker requesting new task
      if(!_mOwnerships.allMapped()) { // there still are tasks
        int mappedSubspcID = _mOwnerships.mappingNextTo(workerID);
        if(mappedSubspcID == ERROR_ID) {
          _prc.abort();
          return pRPL::EVAL_FAILED;
        }
        MPI_Send(&mappedSubspcID, 1, MPI_INT, workerID, pRPL::DSTRBT_TAG, _prc.comm());

        if(readOpt == pRPL::CENT_READING ||
           readOpt == pRPL::CENTDEL_READING) { // centralized reading
          const vector<string> &vInLyrNames = trans.getInLyrNames();
          for(int iLyr = 0; iLyr < vInLyrNames.size(); iLyr++) {
            const string &inLyrName = vInLyrNames[iLyr];
            pRPL::Layer *pLyr = getLayerByName(inLyrName.c_str(), true);
            if(pLyr == NULL) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            if(pLyr->isDecomposed()) { // decomposed Layer
              pRPL::SubCellspace *pSubspc = pLyr->subCellspace_glbID(mappedSubspcID, false);
              if(pSubspc == NULL) {
                if(!pLyr->loadSubCellspaceByGDAL(mappedSubspcID)) {
                  _prc.abort();
                  return pRPL::EVAL_FAILED;
                }
              }
              if(!_iTransfSubspaceBegin(_prc.id(), workerID, inLyrName, mappedSubspcID)) {
                _prc.abort();
                return pRPL::EVAL_FAILED;
              }
            } // end -- if(pLyr->isDecomposed())
          } // end -- for(int iLyr = 0; iLyr < vLyrIDs.size(); iLyr++)
        } // end -- if(readOpt == pRPL::CENT_READING || readOpt == pRPL::CENTDEL_READING)
      } // end -- if(!_mOwnerships.allMapped())

      else { // all tasks are assigned
        pRPL::TaskSignal resp = pRPL::QUIT_SGNL;
        MPI_Send(&resp, 1, MPI_INT, workerID, pRPL::DSTRBT_TAG, _prc.comm());
      }
    } // end -- if((TaskSignal)req == NEWTASK_SGNL)

    else if((pRPL::TaskSignal)req == pRPL::SUBMIT_SGNL) { // worker submitting result
      if(writeOpt == pRPL::PARA_WRITING ||
         writeOpt == pRPL::PARADEL_WRITING) { // parallel writing
        int recvSubspcID;
        MPI_Recv(&recvSubspcID, 1, MPI_INT, workerID, pRPL::SUBMIT_TAG, _prc.comm(), &status);
        if(!mergeTmpSubspcByGDAL(trans, recvSubspcID)) {
          _prc.abort();
          return pRPL::EVAL_FAILED;
        }
      }
      else if(writeOpt == pRPL::CENT_WRITING ||
              writeOpt == pRPL::CENTDEL_WRITING) { // centralized writing
        int recvSubspcInfo[2];
        MPI_Recv(recvSubspcInfo, 2, MPI_INT, workerID, pRPL::SUBMIT_TAG, _prc.comm(), &status);
        int lyrID = recvSubspcInfo[0];
        int subspcGlbID = recvSubspcInfo[1];
        string lyrName = getLayerNameByIdx(lyrID);

        if(!_iTransfSubspaceBegin(workerID, _prc.id(), lyrName, subspcGlbID)) {
          _prc.abort();
          return EVAL_FAILED;
        }
      }
      else if(writeOpt == pRPL::PGT_WRITING ||
              writeOpt == pRPL::PGTDEL_WRITING) {
        int recvSubspcID;
        MPI_Recv(&recvSubspcID, 1, MPI_INT, workerID, pRPL::SUBMIT_TAG, _prc.comm(), &status);
      }
    } // end -- if((TaskSignal)req == SUBMIT_SGNL)

    else if((pRPL::TaskSignal)req == pRPL::QUIT_SGNL) { // worker quitting
      //cout << "Process "<< workerID << " quit" << endl;
      _nActiveWorkers--;
    }

    // test pending transfers, only for centralized I/O
    if(!_vBcstSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltBcstIDs, pRPL::BCST_TRANSF, delSentSpc, writeRecvSpc, delRecvSpc);
    }
    if(!_vSendSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltSendIDs, pRPL::SEND_TRANSF, delSentSpc, writeRecvSpc, delRecvSpc);
    }
    if(!_vRecvSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, delSentSpc, writeRecvSpc, delRecvSpc);
    }

  } // end -- while(_nActiveWorkers > 0)

  // complete all pending transfers, only for centralized I/O
  while(!_vBcstSpcReqs.empty() && nCmpltBcsts != MPI_UNDEFINED) {
    nCmpltBcsts = _iTransfSpaceTest(vCmpltBcstIDs, pRPL::BCST_TRANSF, delSentSpc, writeRecvSpc, delRecvSpc);
  }
  while(!_vSendSpcReqs.empty() && nCmpltSends != MPI_UNDEFINED) {
    nCmpltSends = _iTransfSpaceTest(vCmpltSendIDs, pRPL::SEND_TRANSF, delSentSpc, writeRecvSpc, delRecvSpc);
  }
  while(!_vRecvSpcReqs.empty() && nCmpltRecvs != MPI_UNDEFINED) {
    nCmpltRecvs = _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, delSentSpc, writeRecvSpc, delRecvSpc);
  }

  _clearTransfSpaceReqs();

  return pRPL::EVAL_SUCCEEDED;
}

pRPL::EvaluateReturn pRPL::DataManager::
_writerTF(pRPL::Transition &trans,
          pRPL::WritingOption writeOpt) {
  int nCmpltRecvs = 0;
  pRPL::IntVect vCmpltRecvIDs;
  pRPL::TaskSignal req;
  MPI_Status status;

  bool delSentSpc = false;
  bool writeRecvSpc = (writeOpt != pRPL::NO_WRITING) ? true : false;
  bool delRecvSpc = (writeOpt == pRPL::CENTDEL_WRITING) ? true : false;

  while(_nActiveWorkers > 0) {
    MPI_Recv(&req, 1, MPI_INT, MPI_ANY_SOURCE, pRPL::REQUEST_TAG, _prc.comm(), &status);
    int workerID = status.MPI_SOURCE;

    if((pRPL::TaskSignal)req == pRPL::SUBMIT_SGNL) { // worker submitting result
      if(writeOpt == pRPL::PARA_WRITING ||
         writeOpt == pRPL::PARADEL_WRITING) { // parallel writing
        int recvSubspcID;
        MPI_Recv(&recvSubspcID, 1, MPI_INT, workerID, pRPL::SUBMIT_TAG, _prc.comm(), &status);
        if(!mergeTmpSubspcByGDAL(trans, recvSubspcID)) {
          _prc.abort();
          return pRPL::EVAL_FAILED;
        }
      }
      else if(writeOpt == pRPL::CENT_WRITING ||
              writeOpt == pRPL::CENTDEL_WRITING) { // centralized writing
        int recvSubspcInfo[2];
        MPI_Recv(recvSubspcInfo, 2, MPI_INT, workerID, pRPL::SUBMIT_TAG, _prc.comm(), &status);
        int lyrID = recvSubspcInfo[0];
        int subspcGlbID = recvSubspcInfo[1];
        string lyrName = getLayerNameByIdx(lyrID);

        if(!_iTransfSubspaceBegin(workerID, _prc.id(), lyrName, subspcGlbID)) {
          _prc.abort();
          return EVAL_FAILED;
        }
      }
      else if(writeOpt == pRPL::PGT_WRITING ||
              writeOpt == pRPL::PGTDEL_WRITING) {
        int recvSubspcID;
        MPI_Recv(&recvSubspcID, 1, MPI_INT, workerID, pRPL::SUBMIT_TAG, _prc.comm(), &status);
      }
    } // end -- if((TaskSignal)req == SUBMIT_SGNL)

    else if((pRPL::TaskSignal)req == pRPL::QUIT_SGNL) { // worker quitting
      //cout << "Process "<< workerID << " quit" << endl;
      _nActiveWorkers--;
    }

    // test pending transfers, only for centralized I/O
    if(!_vRecvSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, delSentSpc, writeRecvSpc, delRecvSpc);
    }

  } // end -- while(_nActiveWorkers > 0)

  // complete all pending transfers, only for centralized I/O
  while(!_vRecvSpcReqs.empty() && nCmpltRecvs != MPI_UNDEFINED) {
    nCmpltRecvs = _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, delSentSpc, writeRecvSpc, delRecvSpc);
  }

  _clearTransfSpaceReqs();

  return pRPL::EVAL_SUCCEEDED;
}

pRPL::EvaluateReturn pRPL::DataManager::
_workerTF(pRPL::EvaluateType evalType,
          pRPL::Transition &trans,
          pRPL::ReadingOption readOpt,
          pRPL::WritingOption writeOpt,
          bool ifInitNoData,
          const pRPL::LongVect *pvGlbIdxs) {
  if(evalType != pRPL::EVAL_NONE &&
     evalType != pRPL::EVAL_ALL &&
     evalType != pRPL::EVAL_RANDOMLY &&
     evalType != pRPL::EVAL_SELECTED) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid evaluation type (" \
         << evalType << ")" \
         << endl;
    _prc.abort();
    return pRPL::EVAL_FAILED;
  }

  int nCmpltBcsts = 0, nCmpltSends = 0, nCmpltRecvs = 0;
  pRPL::IntVect vCmpltBcstIDs, vCmpltSendIDs, vCmpltRecvIDs;
  int resp = 0;
  pRPL::TaskSignal req = pRPL::NEWTASK_SGNL;
  MPI_Status status;
  int mappedSubspcID = pRPL::ERROR_ID;
  int submit2ID = _prc.hasWriter() ? 1 : 0; // if a writer exists, send to writer; otherwise master
  int readySubspcID = pRPL::ERROR_ID;
  bool delSentSub = (writeOpt == pRPL::CENTDEL_WRITING) ? true : false;

  while((pRPL::TaskSignal)resp != pRPL::QUIT_SGNL) { // keep running until receiving a Quitting instruction from master
    // if there is less than or equal to one SubCellspace to be evaluated in the queue
    if(_mOwnerships.subspcIDsOnPrc(_prc.id())->size() - _vEvaledSubspcIDs.size() <= 1) {
      req = pRPL::NEWTASK_SGNL;
      MPI_Send(&req, 1, MPI_INT, 0, pRPL::REQUEST_TAG, _prc.comm());
      MPI_Recv(&resp, 1, MPI_INT, 0, pRPL::DSTRBT_TAG, _prc.comm(), &status);
      if((pRPL::TaskSignal)resp != pRPL::QUIT_SGNL) { // received SubCellspace ID from master
        mappedSubspcID = resp;
        _mOwnerships.mappingTo(mappedSubspcID, _prc.id());
        if(!_workerReadSubspcs(trans, mappedSubspcID, readOpt)) {
          _prc.abort();
          return pRPL::EVAL_FAILED;
        }
      } // end -- if((pRPL::TaskSignal)resp != pRPL::QUIT_SGNL)
    } // end -- if(_mOwnerships.subspcIDsOnPrc(_prc.id())->size() - _vEvaledSubspcIDs.size() <= 1)

    // test pending transfers, only for centralized I/O
    if(!_vBcstSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltBcstIDs, pRPL::BCST_TRANSF, delSentSub, false, false);
    }
    if(!_vSendSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltSendIDs, pRPL::SEND_TRANSF, delSentSub, false, false);
    }
    if(!_vRecvSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, delSentSub, false, false);
    }

    // PROCESS A SUBSPACES THAT IS READY
    readySubspcID = _checkTransInDataReady(trans, readOpt);
    if(readySubspcID != pRPL::ERROR_ID) {
      // when task farming, always evaluate the working rectangle.
      // the edge-first method won't be efficient, 'cuz all SubCellspaces may not
      // have been mapped yet.
      //cout << _prc.id() << ": processing SubCellspace " << readySubspcID << endl;
      switch(evalType) {
        case EVAL_NONE:
          if(_evalNone(trans, readySubspcID, ifInitNoData) == pRPL::EVAL_FAILED) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          break;
        case EVAL_ALL:
          if(_evalAll(trans, readySubspcID, pRPL::EVAL_WORKBR, ifInitNoData) == pRPL::EVAL_FAILED) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          break;
        case EVAL_RANDOMLY:
          if(_evalRandomly(trans, readySubspcID, ifInitNoData) == pRPL::EVAL_FAILED) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          break;
        case EVAL_SELECTED:
          if(pvGlbIdxs == NULL) {
            cerr << __FILE__ << " function:" << __FUNCTION__ \
                 << " Error: NULL vector of cell indices to evaluate" \
                 << endl;
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          if(_evalSelected(trans, *pvGlbIdxs, readySubspcID, pRPL::EVAL_WORKBR, ifInitNoData) == pRPL::EVAL_FAILED) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          break;
        default:
          break;
      }
      _vEvaledSubspcIDs.push_back(readySubspcID);

      if(!_workerWriteSubspcs(trans, readySubspcID, writeOpt)) {
        _prc.abort();
        return pRPL::EVAL_FAILED;
      }
    } // end -- if(readySubspcID != ERROR_ID)

  } // end -- while((TaskSignal)resp != QUIT_SGNL)


  // COMPLETE ALL ASSIGNED SUBCELLSPACES
  const pRPL::IntVect *pSubspcIDs = _mOwnerships.subspcIDsOnPrc(_prc.id());
  while(pSubspcIDs != NULL &&
      pSubspcIDs->size() - _vEvaledSubspcIDs.size() > 0) {
    // test pending transfers, only for centralized I/O
    if(!_vBcstSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltBcstIDs, pRPL::BCST_TRANSF, delSentSub, false, false);
    }
    if(!_vSendSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltSendIDs, pRPL::SEND_TRANSF, delSentSub, false, false);
    }
    if(!_vRecvSpcReqs.empty()) {
      _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, delSentSub, false, false);
    }

    // PROCESS A SUBSPACES THAT IS READY
    readySubspcID = _checkTransInDataReady(trans, readOpt);
    if(readySubspcID != pRPL::ERROR_ID) {
      //cout << _prc.id() << ": processing SubCellspace " << readySubspcID << endl;
      switch(evalType) {
        case EVAL_NONE:
          if(_evalNone(trans, readySubspcID, ifInitNoData) == pRPL::EVAL_FAILED) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          break;
        case EVAL_ALL:
          if(_evalAll(trans, readySubspcID, pRPL::EVAL_WORKBR, ifInitNoData) == pRPL::EVAL_FAILED) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          break;
        case EVAL_RANDOMLY:
          if(_evalRandomly(trans, readySubspcID, ifInitNoData) == pRPL::EVAL_FAILED) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          break;
        case EVAL_SELECTED:
          if(pvGlbIdxs == NULL) {
            cerr << __FILE__ << " function:" << __FUNCTION__ \
                 << " Error: NULL vector of cell indices to evaluate" \
                 << endl;
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          if(_evalSelected(trans, *pvGlbIdxs, readySubspcID, pRPL::EVAL_WORKBR, ifInitNoData) == pRPL::EVAL_FAILED) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
          break;
        default:
          break;
      }
      _vEvaledSubspcIDs.push_back(readySubspcID);

      if(!_workerWriteSubspcs(trans, readySubspcID, writeOpt)) {
        _prc.abort();
        return pRPL::EVAL_FAILED;
      }
    } // end -- if(readySubspcID != ERROR_ID)
  }  // end -- while(_mOwnerships.subspcIDsOnPrc(_prc.id())->size() - _vEvaledSubspcIDs.size() > 0)

  // complete all pending transfers, only for centralized I/O
  while(!_vBcstSpcReqs.empty() && nCmpltBcsts != MPI_UNDEFINED) {
    nCmpltBcsts = _iTransfSpaceTest(vCmpltBcstIDs, pRPL::BCST_TRANSF, delSentSub, false, false);
  }
  while(!_vSendSpcReqs.empty() && nCmpltSends != MPI_UNDEFINED) {
    nCmpltSends = _iTransfSpaceTest(vCmpltSendIDs, pRPL::SEND_TRANSF, delSentSub, false, false);
  }
  while(!_vRecvSpcReqs.empty() && nCmpltRecvs != MPI_UNDEFINED) {
    nCmpltRecvs = _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, delSentSub, false, false);
  }

  req = pRPL::QUIT_SGNL;
  MPI_Send(&req, 1, MPI_INT, 0, pRPL::REQUEST_TAG, _prc.comm()); // send a Quitting signal to master
  if(_prc.hasWriter()) {
    MPI_Send(&req, 1, MPI_INT, 1, pRPL::REQUEST_TAG, _prc.comm()); // send a Quitting signal to writer
  }

  _clearTransfSpaceReqs();

  return pRPL::EVAL_SUCCEEDED;
}

pRPL::EvaluateReturn pRPL::DataManager::
_masterST(pRPL::Transition &trans,
          pRPL::WritingOption writeOpt) {
  pRPL::TaskSignal req;
  MPI_Status status;

  if(!_prc.hasWriter()) {
    if(writeOpt == pRPL::PARA_WRITING ||
       writeOpt == pRPL::PARADEL_WRITING) { // parallel writing
      while(_nActiveWorkers > 0) {
        MPI_Recv(&req, 1, MPI_INT, MPI_ANY_SOURCE, pRPL::REQUEST_TAG, _prc.comm(), &status);
        int workerID = status.MPI_SOURCE;
        if((pRPL::TaskSignal)req == pRPL::QUIT_SGNL) { // worker quitting
          //cout << "Process "<< workerID << " quit" << endl;
          _nActiveWorkers--;
        }
      } // end -- while(_nActiveWorkers > 0)
      if(!mergeAllTmpSubspcsByGDAL(trans)) {
        _prc.abort();
        return pRPL::EVAL_FAILED;;
      }
    } // end -- parallel writing

    else if(writeOpt == pRPL::CENT_WRITING ||
            writeOpt == pRPL::CENTDEL_WRITING) { // centralized writing
      bool delRecvSpc = (writeOpt == pRPL::CENTDEL_WRITING) ? true : false;
      int nCmpltRecvs = 0;
      pRPL::IntVect vCmpltRecvIDs;

      while(_nActiveWorkers > 0) {
        MPI_Recv(&req, 1, MPI_INT, MPI_ANY_SOURCE, pRPL::REQUEST_TAG, _prc.comm(), &status);
        int workerID = status.MPI_SOURCE;
        if((pRPL::TaskSignal)req == pRPL::SUBMIT_SGNL) { // worker submitting result
          int recvSubspcInfo[2];
          MPI_Recv(recvSubspcInfo, 2, MPI_INT, workerID, pRPL::SUBMIT_TAG, _prc.comm(), &status);
          int lyrID = recvSubspcInfo[0];
          int subspcGlbID = recvSubspcInfo[1];
          string lyrName = getLayerNameByIdx(lyrID);

          if(!_iTransfSubspaceBegin(workerID, _prc.id(), lyrName, subspcGlbID)) {
            _prc.abort();
            return EVAL_FAILED;
          }
        } // end -- if((TaskSignal)req == SUBMIT_SGNL)

        else if((pRPL::TaskSignal)req == pRPL::QUIT_SGNL) { // worker quitting
          //cout << "Process "<< workerID << " quit" << endl;
          _nActiveWorkers--;
        }

        if(!_vRecvSpcReqs.empty()) {
          _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, true, true, delRecvSpc);
        }
      } // end -- while(_nActiveWorkers > 0)

      // write local SubCellspaces
      const IntVect *pvSubspcIDs = _mOwnerships.subspcIDsOnPrc(_prc.id());
      if(pvSubspcIDs != NULL) {
        for(int iSubspc = 0; iSubspc < pvSubspcIDs->size(); iSubspc++) {
          int subspcGlbID = pvSubspcIDs->at(iSubspc);
          const vector<string> &vOutLyrNames = trans.getOutLyrNames();
          for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
            const string &outLyrName = vOutLyrNames[iOutLyr];
            pRPL::Layer *pOutLyr = getLayerByName(outLyrName.c_str(), true);
            if(pOutLyr == NULL) {
              _prc.abort();
              return EVAL_FAILED;
            }
            if(!pOutLyr->writeSubCellspaceByGDAL(subspcGlbID)) {
              _prc.abort();
              return EVAL_FAILED;
            }
          } // end -- for(iOutLyr)
        } // end -- for(iSubspc)
      }

      while(!_vRecvSpcReqs.empty() && nCmpltRecvs != MPI_UNDEFINED) {
        nCmpltRecvs = _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, true, true, delRecvSpc);
      }

    } // end -- centralized writing

    else if(writeOpt == pRPL::PGT_WRITING ||
            writeOpt == pRPL::PGTDEL_WRITING) { // parallel writing
      while(_nActiveWorkers > 0) {
        MPI_Recv(&req, 1, MPI_INT, MPI_ANY_SOURCE, pRPL::REQUEST_TAG, _prc.comm(), &status);
        int workerID = status.MPI_SOURCE;
        if((pRPL::TaskSignal)req == pRPL::QUIT_SGNL) { // worker quitting
          //cout << "Process "<< workerID << " quit" << endl;
          _nActiveWorkers--;
        }
      } // end -- while(_nActiveWorkers > 0)
    } // end -- parallel writing
  }  // end -- if(!_prc.hasWriter())

  _clearTransfSpaceReqs();

  return pRPL::EVAL_SUCCEEDED;
}

pRPL::EvaluateReturn pRPL::DataManager::
_writerST(pRPL::Transition &trans,
          pRPL::WritingOption writeOpt) {
  return _writerTF(trans, writeOpt);
}

pRPL::EvaluateReturn pRPL::DataManager::
_workerST(pRPL::EvaluateType evalType,
          pRPL::Transition &trans,
          pRPL::WritingOption writeOpt,
          bool ifInitNoData,
          const pRPL::LongVect *pvGlbIdxs) {
  if(evalType != pRPL::EVAL_NONE &&
     evalType != pRPL::EVAL_ALL &&
     evalType != pRPL::EVAL_RANDOMLY &&
     evalType != pRPL::EVAL_SELECTED) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid evaluation type (" \
         << evalType << ")" \
         << endl;
    _prc.abort();
    return pRPL::EVAL_FAILED;
  }

  _clearTransfSpaceReqs();

  const IntVect *pvSubspcIDs = _mOwnerships.subspcIDsOnPrc(_prc.id());
  if(pvSubspcIDs != NULL) {
    for(int iSubspc = 0; iSubspc < pvSubspcIDs->size(); iSubspc++) {
      int subspcGlbID = pvSubspcIDs->at(iSubspc);
      if(!trans.needExchange() ||
         (trans.needExchange() && !trans.edgesFirst())) { // no exchange needed
        //cout << _prc.id() << ": processing SubCellspace " << subspcGlbID << endl;
        switch(evalType) {
          case EVAL_NONE:
            if(_evalNone(trans, subspcGlbID, ifInitNoData) == pRPL::EVAL_FAILED) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            break;
          case EVAL_ALL:
            if(_evalAll(trans, subspcGlbID, pRPL::EVAL_WORKBR, ifInitNoData) == pRPL::EVAL_FAILED) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            break;
          case EVAL_RANDOMLY:
            if(_evalRandomly(trans, subspcGlbID, ifInitNoData) == pRPL::EVAL_FAILED) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            break;
          case EVAL_SELECTED:
            if(pvGlbIdxs == NULL) {
              cerr << __FILE__ << " function:" << __FUNCTION__ \
                   << " Error: NULL vector of cell indices to evaluate" \
                   << endl;
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            if(_evalSelected(trans, *pvGlbIdxs, subspcGlbID, pRPL::EVAL_WORKBR, ifInitNoData) == pRPL::EVAL_FAILED) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            break;
          default:
            break;
        }
        _vEvaledSubspcIDs.push_back(subspcGlbID);

        if(!_dymWriteSubspcs(trans, subspcGlbID, writeOpt)) {
          _prc.abort();
          return pRPL::EVAL_FAILED;
        }
      } // end -- if(!trans.needExchange())

      else { // exchange needed, process the edges first
        //cout << _prc.id() << ": processing SubCellspace edge " << subspcGlbID << endl;
        switch(evalType) {
          case EVAL_NONE:
            if(_evalNone(trans, subspcGlbID, ifInitNoData) == pRPL::EVAL_FAILED) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            break;
          case EVAL_ALL:
            if(_evalAll(trans, subspcGlbID, pRPL::EVAL_EDGES, ifInitNoData) == pRPL::EVAL_FAILED) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            break;
          case EVAL_RANDOMLY:
            if(_evalRandomly(trans, subspcGlbID, ifInitNoData) == pRPL::EVAL_FAILED) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            break;
          case EVAL_SELECTED:
            if(pvGlbIdxs == NULL) {
              cerr << __FILE__ << " function:" << __FUNCTION__ \
                  << " Error: NULL vector of cell indices to evaluate" \
                  << endl;
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            if(_evalSelected(trans, *pvGlbIdxs, subspcGlbID, pRPL::EVAL_EDGES, ifInitNoData) == pRPL::EVAL_FAILED) {
              _prc.abort();
              return pRPL::EVAL_FAILED;
            }
            break;
          default:
            break;
        }
      } // end -- exchange needed
    } // end -- for(iSubspc)

    if(trans.needExchange()) { // exchange, while processing the interiors
      if(!_makeCellStream(trans) || !_iExchangeBegin()) {
        _prc.abort();
        return pRPL::EVAL_FAILED;
      }

      if(trans.edgesFirst()) {
        for(int iSubspc = 0; iSubspc < pvSubspcIDs->size(); iSubspc++) {
          int subspcGlbID = pvSubspcIDs->at(iSubspc);
          //cout << _prc.id() << ": processing SubCellspace interior " << subspcGlbID << endl;
          switch(evalType) {
            case EVAL_NONE:
              break;
            case EVAL_ALL:
              if(_evalAll(trans, subspcGlbID, pRPL::EVAL_INTERIOR, ifInitNoData) == pRPL::EVAL_FAILED) {
                _prc.abort();
                return pRPL::EVAL_FAILED;
              }
              break;
            case EVAL_RANDOMLY:
              break;
            case EVAL_SELECTED:
              if(pvGlbIdxs == NULL) {
                cerr << __FILE__ << " function:" << __FUNCTION__ \
                    << " Error: NULL vector of cell indices to evaluate" \
                    << endl;
                _prc.abort();
                return pRPL::EVAL_FAILED;
              }
              if(_evalSelected(trans, *pvGlbIdxs, subspcGlbID, pRPL::EVAL_INTERIOR, ifInitNoData) == pRPL::EVAL_FAILED) {
                _prc.abort();
                return pRPL::EVAL_FAILED;
              }
              break;
            default:
              break;
          }
          _vEvaledSubspcIDs.push_back(subspcGlbID);
          if(!_dymWriteSubspcs(trans, subspcGlbID, writeOpt)) {
            _prc.abort();
            return pRPL::EVAL_FAILED;
          }
        } // end -- for(iSubspc)

        if(!_iExchangeEnd()) {
          return pRPL::EVAL_FAILED;
        }
      }
    } // end -- if(trans.needExchange())

    if(!_finalWriteSubspcs(trans, writeOpt)) {
      _prc.abort();
      return pRPL::EVAL_FAILED;
    }

  } // end -- if(pvSubspcIDs != NULL)

  // complete all pending transfers, only for centralized I/O
  int nCmpltBcsts = 0, nCmpltSends = 0, nCmpltRecvs = 0;
  pRPL::IntVect vCmpltBcstIDs, vCmpltSendIDs, vCmpltRecvIDs;
  bool delSentSub = (writeOpt == pRPL::CENTDEL_WRITING) ? true : false;

  while(!_vSendSpcReqs.empty() && nCmpltSends != MPI_UNDEFINED) {
    nCmpltSends = _iTransfSpaceTest(vCmpltSendIDs, pRPL::SEND_TRANSF, delSentSub, false, false);
  }

  // send Quitting signal
  pRPL::TaskSignal req = pRPL::QUIT_SGNL;
  if(_prc.hasWriter()) {
    MPI_Send(&req, 1, MPI_INT, 1, pRPL::REQUEST_TAG, _prc.comm()); // send a Quitting signal to writer
  }

  if(!_prc.isMaster()) {
    if(!_prc.hasWriter() && writeOpt != pRPL::NO_WRITING) {
      MPI_Send(&req, 1, MPI_INT, 0, pRPL::REQUEST_TAG, _prc.comm()); // send a Quitting signal to master
      //cout << _prc.id() << " quitting" << endl;
    }
    _clearTransfSpaceReqs();
  }

  return pRPL::EVAL_SUCCEEDED;
}

bool pRPL::DataManager::
_calcExchangeBRs(pRPL::Transition &trans) {
  const vector<string> &vOutLyrNames = trans.getOutLyrNames();
  for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
    pRPL::Layer *pOutLyr = getLayerByName(vOutLyrNames[iOutLyr].c_str(), true);
    if(pOutLyr == NULL) {
      return false;
    }
    if(!pOutLyr->calcAllBRs(trans.onlyUpdtCtrCell(), getNbrhdByName(trans.getNbrhdName().c_str()))) {
      return false;
    }
  }
  return true;
}

bool pRPL::DataManager::
_calcExchangeRoutes(pRPL::Transition &trans) {
  _mSendRoutes.clear();
  _mRecvRoutes.clear();

  const vector<string> &vOutLyrNames = trans.getOutLyrNames();
  const pRPL::IntVect *pvAssignments = _mOwnerships.subspcIDsOnPrc(_prc.id());

  if(trans.needExchange() && !vOutLyrNames.empty() &&
     pvAssignments != NULL && !pvAssignments->empty()) {
    pRPL::Layer *pOutLyr = getLayerByName(vOutLyrNames[0].c_str(), true); // only one output Layer is needed for calculating the exchange map
    if(pOutLyr == NULL) {
      return false;
    }
    for(int iAssignment = 0; iAssignment < pvAssignments->size(); iAssignment++) {
      int subspcID = pvAssignments->at(iAssignment);
      const pRPL::SubCellspaceInfo *pSubspcInfo = pOutLyr->subCellspaceInfo_glbID(subspcID, true);
      if(pSubspcInfo == NULL) {
        return false;
      }

      for(int iDir = 0; iDir < pSubspcInfo->nNbrDirs(); iDir++) {
        const pRPL::IntVect &vNbrSpcIDs = pSubspcInfo->nbrSubSpcIDs(iDir);
        for(int iNbr = 0; iNbr < vNbrSpcIDs.size(); iNbr++) {
          int nbrSubspcID = vNbrSpcIDs[iNbr];
          int nbrPrc = _mOwnerships.findPrcBySubspcID(nbrSubspcID);

          const pRPL::CoordBR *pSendBR = pSubspcInfo->sendBR(iDir, iNbr);
          if(pSendBR != NULL && nbrPrc >= 0) {
            if(!_mSendRoutes.add(subspcID, iDir, iNbr, nbrPrc)) {
              return false;
            }
          }

          const pRPL::SubCellspaceInfo *pNbrSubspcInfo = pOutLyr->subCellspaceInfo_glbID(nbrSubspcID, true);
          if(pNbrSubspcInfo == NULL) {
            return false;
          }
          int iOppDir = pNbrSubspcInfo->oppositeDir(iDir);
          const pRPL::IntVect &vOppNbrSpcIDs = pNbrSubspcInfo->nbrSubSpcIDs(iOppDir);
          pRPL::IntVect::const_iterator itrOppNbrSpcID = std::find(vOppNbrSpcIDs.begin(), vOppNbrSpcIDs.end(), subspcID);
          if(itrOppNbrSpcID != vOppNbrSpcIDs.end()) {
            int iOppNbr = itrOppNbrSpcID - vOppNbrSpcIDs.begin();
            const pRPL::CoordBR *pOppSendBR = pNbrSubspcInfo->sendBR(iOppDir, iOppNbr);
            if(pOppSendBR != NULL && nbrPrc >= 0) {
              if(!_mRecvRoutes.add(subspcID, iDir, iNbr, nbrPrc)) {
                return false;
              }
            }
          }
        } // end -- for(iNbr)
      } // end -- for(iDir)
    } // end -- for(iAssignment)
  }

  return true;
}

bool pRPL::DataManager::
_loadCellStream(pRPL::CellStream &cells2load) {
  long cellGlbIdx;
  void *aCellVal;
  const vector<pair<unsigned int, unsigned int> >& vLyrInfos = cells2load.getLayerInfos();
  for(int iLyr = 0; iLyr < vLyrInfos.size(); iLyr++) {
    unsigned int lyrID = vLyrInfos[iLyr].first;
    unsigned int dataSize = vLyrInfos[iLyr].second;
    int nCells = cells2load.getNumCellsOnLayer(lyrID);
    void* aCells = cells2load.getCellsOnLayer(lyrID);

    pRPL::Layer *pLyr = getLayerByIdx(lyrID, true);
    if(pLyr == NULL) {
      return false;
    }
    if(!pLyr->loadCellStream(aCells, nCells)) {
      return false;
    }
  } // end --- for(iLyr)
  return true;
}

bool pRPL::DataManager::
_makeCellStream(pRPL::Transition &trans) {
  _vSendCellReqs.clear();
  _vRecvCellReqs.clear();
  _mSendCells.clear();
  _mRecvCells.clear();

  const vector<string> &vOutLyrNames = trans.getOutLyrNames();
  if(vOutLyrNames.empty()) {
    return true;
  }

  pRPL::ExchangeMap::iterator itrSendRoute = _mSendRoutes.begin();
  while(itrSendRoute != _mSendRoutes.end()) {
    int toPrc = itrSendRoute->first;
    vector<pRPL::ExchangeNode> &vSendNodes = itrSendRoute->second;
    for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
      pRPL::Layer *pOutLyr = getLayerByName(vOutLyrNames[iOutLyr].c_str(), true);
      if(pOutLyr == NULL) {
        return false;
      }
      int lyrID = getLayerIdxByName(vOutLyrNames[iOutLyr].c_str(), true);
      unsigned int lyrDataSize = pOutLyr->dataSize();
      if(!_mSendCells[toPrc].addLayer(lyrID, lyrDataSize)) {
        return false;
      }

      for(int iNode = 0; iNode < vSendNodes.size(); iNode++) {
        int subspcID = vSendNodes[iNode].subspcGlbID;
        int iDir = vSendNodes[iNode].iDirection;
        int iNbr = vSendNodes[iNode].iNeighbor;
        pRPL::SubCellspace *pSubspc = pOutLyr->subCellspace_glbID(subspcID, true);
        if(pSubspc == NULL) {
          return false;
        }
        const pRPL::CoordBR *pSendBR = pSubspc->subInfo()->sendBR(iDir, iNbr);
        const pRPL::LongVect &vUpdtIdxs = pSubspc->getUpdatedIdxs();
        for(int iIdx = 0; iIdx < vUpdtIdxs.size(); iIdx++) {
          long lclIdx = vUpdtIdxs[iIdx];
          if(pSendBR->ifContain(lclIdx, pSubspc->info()->dims())) {
            long glbIdx = pSubspc->subInfo()->lclIdx2glbIdx(lclIdx);
            void *aCellVal = pSubspc->at(lclIdx);
            if(!_mSendCells[toPrc].addCell(glbIdx, aCellVal)) {
              return false;
            }
            //cout << toPrc << " " << pOutLyr->name() << " [" << glbIdx << ", " << pSubspc->atAs<unsigned short>(lclIdx) << "]" << endl;
          }
        } // end -- for(iIdx)
      } // end -- for(iNode)
    } // end -- for(iOutLyr)

    if(toPrc != _prc.id()) { // if toPrc is not myself
      _vSendCellReqs.resize(_vSendCellReqs.size() + 1);
      MPI_Isend(&(_mSendCells[toPrc].getCellCounts().at(0)),
                _mSendCells[toPrc].getNumLayers(),
                MPI_INT, toPrc, _prc.id(), _prc.comm(),
                &(_vSendCellReqs.back()));
    }

    itrSendRoute++;
  } // end -- while(itrSendRoute != _mSendRoutes.end())

  // clear all update tracks
  for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
    pRPL::Layer *pOutLyr = getLayerByName(vOutLyrNames[iOutLyr].c_str(), true);
    if(pOutLyr == NULL) {
      return false;
    }
    pOutLyr->clearUpdateTracks();
  }

  pRPL::ExchangeMap::iterator itrRecvRoute = _mRecvRoutes.begin();
  while(itrRecvRoute != _mRecvRoutes.end()) {
    int fromPrc = itrRecvRoute->first;
    vector<pRPL::ExchangeNode> &vRecvNodes = itrRecvRoute->second;
    for(int iOutLyr = 0; iOutLyr < vOutLyrNames.size(); iOutLyr++) {
      pRPL::Layer *pOutLyr = getLayerByName(vOutLyrNames[iOutLyr].c_str(), true);
      if(pOutLyr == NULL) {
        return false;
      }
      int lyrID = getLayerIdxByName(vOutLyrNames[iOutLyr].c_str(), true);
      unsigned int lyrDataSize = pOutLyr->dataSize();
      if(!_mRecvCells[fromPrc].addLayer(lyrID, lyrDataSize)) {
        return false;
      }
    } // end -- for(iOutLyr)

    if(fromPrc != _prc.id()) {
      _vRecvCellReqs.resize(_vRecvCellReqs.size() + 1);
      MPI_Irecv(&(_mRecvCells[fromPrc].getCellCounts().at(0)),
                _mRecvCells[fromPrc].getNumLayers(),
                MPI_INT, fromPrc, fromPrc, _prc.comm(),
                &(_vRecvCellReqs.back()));
    }

    itrRecvRoute++;
  } // end -- while(itrRecvRoute != _mRecvRoutes.end())

  return true;
}

bool pRPL::DataManager::
_iExchangeBegin() {
  vector<MPI_Status> vSendStats(_vSendCellReqs.size());
  MPI_Waitall(_vSendCellReqs.size(), &(_vSendCellReqs[0]), &(vSendStats[0]));
  _vSendCellReqs.clear();

  pRPL::ExchangeMap::iterator itrSendRoute = _mSendRoutes.begin();
  while(itrSendRoute != _mSendRoutes.end()) {
    int toPrc = itrSendRoute->first;
    if(toPrc != _prc.id()) {
      _vSendCellReqs.resize(_vSendCellReqs.size() + 1);
      MPI_Isend(_mSendCells[toPrc].getStream(),
                _mSendCells[toPrc].size(),
                MPI_CHAR, toPrc, _prc.id(), _prc.comm(),
                &(_vSendCellReqs.back()));
    }

    itrSendRoute++;
  } // end -- while(itrSendRoute != _mSendRoutes.end())

  vector<MPI_Status> vRecvStats(_vRecvCellReqs.size());
  MPI_Waitall(_vRecvCellReqs.size(), &(_vRecvCellReqs[0]), &(vRecvStats[0]));
  _vRecvCellReqs.clear();

  pRPL::ExchangeMap::iterator itrRecvRoute = _mRecvRoutes.begin();
  while(itrRecvRoute != _mRecvRoutes.end()) {
    int fromPrc = itrRecvRoute->first;
    if(fromPrc != _prc.id()) {
      if(!_mRecvCells[fromPrc].resize()) {
        return false;
      }
      _vRecvCellReqs.resize(_vRecvCellReqs.size() + 1);
      MPI_Irecv(_mRecvCells[fromPrc].getStream(),
                _mRecvCells[fromPrc].size(),
                MPI_CHAR, fromPrc, fromPrc, _prc.comm(),
                &(_vRecvCellReqs.back()));
    }
    else {
      _mRecvCells[_prc.id()] = _mSendCells[_prc.id()];
    }

    itrRecvRoute++;
  } // end --- while(itrRecvRoute != _mRecvRoutes.end())

  return true;
}

bool pRPL::DataManager::
_iExchangeEnd() {
  vector<MPI_Status> vSendStats(_vSendCellReqs.size());
  MPI_Waitall(_vSendCellReqs.size(), &(_vSendCellReqs[0]), &(vSendStats[0]));
  _vSendCellReqs.clear();
  _mSendCells.clear();

  vector<MPI_Status> vRecvStats(_vRecvCellReqs.size());
  MPI_Waitall(_vRecvCellReqs.size(), &(_vRecvCellReqs[0]), &(vRecvStats[0]));
  _vRecvCellReqs.clear();
  map<int, pRPL::CellStream>::iterator itrRecvCells = _mRecvCells.begin();
  while(itrRecvCells != _mRecvCells.end()) {
    pRPL::CellStream &recvCells = itrRecvCells->second;
    if(!_loadCellStream(recvCells)) {
      return false;
    }
    itrRecvCells++;
  }
  _mRecvCells.clear();

  return true;
}

/****************************************************
*                 Public Methods                    *
****************************************************/
pRPL::DataManager::
DataManager()
  :_nActiveWorkers(0) {
  GDALAllRegister();
}

pRPL::DataManager::
~DataManager() {
  clearLayers();
  clearNbrhds();
}

bool pRPL::DataManager::
initMPI(MPI_Comm comm,
        bool hasWriter) {
  return _prc.set(comm, hasWriter);
}

void pRPL::DataManager::
finalizeMPI() {
  _prc.finalize();
}

pRPL::Process& pRPL::DataManager::
mpiPrc() {
  return _prc;
}

int pRPL::DataManager::
nLayers() const {
  return _lLayers.size();
}

pRPL::Layer* pRPL::DataManager::
addLayer(const char *aLyrName) {
  pRPL::Layer *pLyr = getLayerByName(aLyrName, false);
  if(pLyr != NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__
         << " Error: Layer name (" << aLyrName << ") already exists in DataManager. " \
         << "Failed to add new Layer."<< endl;
  }
  else {
    
    _lLayers.push_back(pRPL::Layer(aLyrName));
    pLyr = &(_lLayers.back());
   
  }
  return pLyr;
}


/*---------------------------GDAL----------------------------------------------------*/
pRPL::Layer* pRPL::DataManager::
addLayerByGDAL(const char *aLyrName,
               const char *aGdalFileName,
               int iBand,
               bool pReading) {
  pRPL::Layer *pLyr = addLayer(aLyrName);

  if(pReading) { // parallel reading: all processes read the file in parallel
    if(!pLyr->openGdalDS(aGdalFileName, GA_ReadOnly)) {
      _prc.abort();
      return pLyr;
    }
    pLyr->dsBand(iBand);

    if(!pLyr->initCellspaceInfoByGDAL()) {
      _prc.abort();
      return pLyr;
    }
  }
  else { // master process reads the file and distribute info to other processes
    vector<char> vBuf;
    int bufSize;

    if(_prc.isMaster()) {
      if(!pLyr->openGdalDS(aGdalFileName, GA_ReadOnly)) {
        _prc.abort();
        return pLyr;
      }
      pLyr->dsBand(iBand);

      if(!pLyr->initCellspaceInfoByGDAL()) {
        _prc.abort();
        return pLyr;
      }

      pLyr->cellspaceInfo()->add2Buf(vBuf);
      bufSize = vBuf.size();
    }
    else {
      pLyr->initCellspaceInfo();
      pLyr->dsBand(iBand); // worker will know this Layer has a GDAL dataset
    }

    if(MPI_Bcast(&bufSize, 1, MPI_INT, 0, _prc.comm()) != MPI_SUCCESS) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: failed to broadcast CellspaceInfo buffer size" \
          << endl;
      _prc.abort();
      return pLyr;
    }

    if(!_prc.isMaster()) {
      vBuf.resize(bufSize);
    }

    if(MPI_Bcast(&(vBuf[0]), bufSize, MPI_CHAR, 0, _prc.comm()) != MPI_SUCCESS) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: failed to broadcast CellspaceInfo" \
          << endl;
      _prc.abort();
      return pLyr;
    }

    if(!_prc.isMaster()) {
      if(!pLyr->cellspaceInfo()->initFromBuf(vBuf)) {
        _prc.abort();
        return pLyr;
      }
    }
  }
  return pLyr;
}

bool pRPL::DataManager::
createLayerGDAL(const char *aLyrName,
                const char *aGdalFileName,
                const char *aGdalFormat,
                char **aGdalOptions) {
  pRPL::Layer *pLyr = getLayerByName(aLyrName, true);
  if(pLyr == NULL) {
    _prc.abort();
    return false;
  }

  if(_prc.hasWriter()) { // if a writer Process exists
    if(_prc.isWriter()) {
      if(!pLyr->createGdalDS(aGdalFileName, aGdalFormat, aGdalOptions)) {
        _prc.abort();
        return false;
      }
    }
    _prc.sync();
    if(!_prc.isWriter()) {
      pLyr->openGdalDS(aGdalFileName, GA_ReadOnly); // worker processes cannot write to the dataset
      pLyr->dsBand(1); // all Processes will know this layer has a GDAL dataset
    }
  }
  else { // if no writer Process exists
    if(_prc.isMaster()) {
      if(!pLyr->createGdalDS(aGdalFileName, aGdalFormat, aGdalOptions)) {
        _prc.abort();
        return false;
      }
    }
    _prc.sync();
    if(!_prc.isMaster()) {
      pLyr->openGdalDS(aGdalFileName, GA_ReadOnly); // worker processes cannot write to the dataset
      pLyr->dsBand(1); // all Processes will know this layer has a GDAL dataset
    }
  }

  return true;
}

void pRPL::DataManager::
closeGDAL() {
  if(!_lLayers.empty()) {
    list<pRPL::Layer>::iterator itrLyr = _lLayers.begin();
    while(itrLyr != _lLayers.end()) {
      itrLyr->closeGdalDS();
      itrLyr++;
    }
  }
}

/*---------------------------hsj-------------PGIOL--------------------*/
pRPL::Layer* pRPL::DataManager::
addLayerByPGTIOL(const char *aLyrName,
                 const char *aPgtiolFileName,
                 int iBand,
                 bool pReading) {
  pRPL::Layer *pLyr = addLayer(aLyrName);

  if(pReading) { // parallel reading: all processes read the file in parallel
    int prc = _prc.hasWriter()?1:0;
    if(!pLyr->openPgtiolDS(aPgtiolFileName, prc, PG_ReadOnly)) {
      _prc.abort();
      return pLyr;
    }
    pLyr->dsBand(iBand);

    if(!pLyr->initCellspaceInfoByPGTIOL()) {
      _prc.abort();
      return pLyr;
    }
  }
  else { // master process reads the file and distribute info to other processes
    vector<char> vBuf;
    int bufSize;

    if(_prc.isMaster()) {
      if(!pLyr->openPgtiolDS(aPgtiolFileName,_prc.id())) {
        _prc.abort();
        return pLyr;
      }
      pLyr->dsBand(iBand);

      if(!pLyr->initCellspaceInfoByPGTIOL()) {
        _prc.abort();
        return pLyr;
      }

      pLyr->cellspaceInfo()->add2Buf(vBuf);
      bufSize = vBuf.size();
    }
    else {
      pLyr->initCellspaceInfo();
      pLyr->dsBand(iBand); // worker will know this Layer has a GDAL dataset
    }

    if(MPI_Bcast(&bufSize, 1, MPI_INT, 0, _prc.comm()) != MPI_SUCCESS) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: failed to broadcast CellspaceInfo buffer size" \
          << endl;
      _prc.abort();
      return pLyr;
    }

    if(!_prc.isMaster()) {
      vBuf.resize(bufSize);
    }

    if(MPI_Bcast(&(vBuf[0]), bufSize, MPI_CHAR, 0, _prc.comm()) != MPI_SUCCESS) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: failed to broadcast CellspaceInfo" \
          << endl;
      _prc.abort();
      return pLyr;
    }

    if(!_prc.isMaster()) {
      if(!pLyr->cellspaceInfo()->initFromBuf(vBuf)) {
        _prc.abort();
        return pLyr;
      }
    }
  }
  return pLyr;
}

bool pRPL::DataManager::
createLayerPGTIOL(const char *aLyrName,
                  const char *aPgtiolFileName,
                  char **aPgtiolOptions) {
  pRPL::Layer *pLyr = getLayerByName(aLyrName, true);
  if(pLyr == NULL) {
    _prc.abort();
    return false;
  } 
  if(_prc.hasWriter()) { // if a writer Process exists
      if(_prc.isWriter()) {
      if(!pLyr->createPgtiolDS(aPgtiolFileName,1,aPgtiolOptions)) {
        _prc.abort();
        return false;
      }
    }
    _prc.sync();
    if(!_prc.isWriter()) {
      pLyr->openPgtiolDS(aPgtiolFileName,1,PG_Update);
      pLyr->dsBand(1); // all Processes will know this layer has a pGTIOL dataset
    }
  }
  else { // if no writer Process exists
    
    if(_prc.isMaster()) {
      if(!pLyr->createPgtiolDS(aPgtiolFileName,0,aPgtiolOptions)) {
        _prc.abort();
        return false;
      }
    }
	 
    _prc.sync();
    if(!_prc.isMaster()) {
      pLyr->openPgtiolDS(aPgtiolFileName,0,PG_Update);
      pLyr->dsBand(1); // all Processes will know this layer has a pGTIOL dataset
    }
  }
  return true;
}

void pRPL::DataManager::
closePGTIOL() {
  if(!_lLayers.empty()) {
    list<pRPL::Layer>::iterator itrLyr = _lLayers.begin();
    while(itrLyr != _lLayers.end()) {
      itrLyr->closePgtiolDS();
      itrLyr++;
    }
  }
}

void pRPL::DataManager::
closeDatasets() {
  if(!_lLayers.empty()) {
    list<pRPL::Layer>::iterator itrLyr = _lLayers.begin();
    while(itrLyr != _lLayers.end()) {
      itrLyr->closeGdalDS();
      itrLyr->closePgtiolDS();
      itrLyr++;
    }
  }
}

bool pRPL::DataManager::
rmvLayerByIdx(int lyrID,
              bool warning) {
  bool done = true;
  if(lyrID < 0 || lyrID >= _lLayers.size()) {
    done = false;
    if(warning) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: invalid Layer index (" << lyrID \
           << ")" << endl;
    }
  }
  else {
    list<pRPL::Layer>::iterator itrLyr = _lLayers.begin();
    advance(itrLyr, lyrID);
    _lLayers.erase(itrLyr);
  }
  return done;
}

bool pRPL::DataManager::
rmvLayerByName(const char *aLyrName,
               bool warning) {
  bool done = false;
  list<pRPL::Layer>::iterator itrLyr = _lLayers.begin();
  while(itrLyr != _lLayers.end()) {
    if(strcmp(aLyrName, itrLyr->name()) == 0) {
      _lLayers.erase(itrLyr);
      done = true;
      break;
    }
    itrLyr++;
  }
  if(!done && warning) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: invalid Layer name (" << aLyrName \
        << ")" << endl;
  }
  return done;
}

void pRPL::DataManager::
clearLayers() {
  _lLayers.clear();
}

pRPL::Layer* pRPL::DataManager::
getLayerByIdx(int lyrID,
              bool warning) {
  pRPL::Layer *pLayer = NULL;
  if(lyrID >= 0 && lyrID < _lLayers.size()) {
    list<pRPL::Layer>::iterator itrLyr = _lLayers.begin();
    advance(itrLyr, lyrID);
    pLayer = &(*itrLyr);
  }
  else {
    if(warning) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Warning: NO Layer with index (" << lyrID \
           << ") was found" << endl;
    }
  }
  return pLayer;
}

const pRPL::Layer* pRPL::DataManager::
getLayerByIdx(int lyrID,
              bool warning) const {
  const pRPL::Layer *pLayer = NULL;
  if(lyrID >= 0 && lyrID < _lLayers.size()) {
    list<pRPL::Layer>::const_iterator itrLyr = _lLayers.begin();
    advance(itrLyr, lyrID);
    pLayer = &(*itrLyr);
  }
  else {
    if(warning) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Warning: NO Layer with index (" << lyrID \
           << ") was found" << endl;
    }
  }
  return pLayer;
}

int pRPL::DataManager::
getLayerIdxByName(const char *aLyrName,
                  bool warning) const {
  int lyrIdx = ERROR_ID, lyrCount = 0;
  bool found = false;
  list<pRPL::Layer>::const_iterator itrLyr = _lLayers.begin();
  while(itrLyr != _lLayers.end()) {
    if(strcmp(aLyrName, itrLyr->name()) == 0) {
      lyrIdx = lyrCount;
      found = true;
      break;
    }
    lyrCount++;
    itrLyr++;
  }
  if(!found && warning) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Warning: NO Layer with name (" << aLyrName \
         << ") was found" << endl;
  }
  return lyrIdx;
}

pRPL::Layer* pRPL::DataManager::
getLayerByName(const char *aLyrName,
               bool warning) {
  pRPL::Layer* pLyr = NULL;
  list<pRPL::Layer>::iterator itrLyr = _lLayers.begin();
  while(itrLyr != _lLayers.end()) {
    if(strcmp(aLyrName, itrLyr->name()) == 0) {
      pLyr = &(*itrLyr);
      break;
    }
    itrLyr++;
  }
  if(pLyr == NULL && warning) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Warning: NO Layer with name (" << aLyrName \
        << ") was found" << endl;
  }
  return pLyr;
}

const pRPL::Layer* pRPL::DataManager::
getLayerByName(const char *aLyrName,
               bool warning) const {
  const pRPL::Layer* pLyr = NULL;
  list<pRPL::Layer>::const_iterator itrLyr = _lLayers.begin();
  while(itrLyr != _lLayers.end()) {
    if(strcmp(aLyrName, itrLyr->name()) == 0) {
      pLyr = &(*itrLyr);
      break;
    }
    itrLyr++;
  }
  if(pLyr == NULL && warning) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Warning: NO Layer with name (" << aLyrName \
        << ") was found" << endl;
  }
  return pLyr;
}

const char* pRPL::DataManager::
getLayerNameByIdx(int lyrID,
                  bool warning) const {
  const char *aLyrName = NULL;
  const pRPL::Layer *pLyr = getLayerByIdx(lyrID, warning);
  if(pLyr != NULL) {
    aLyrName = pLyr->name();
  }
  return aLyrName;
}

bool pRPL::DataManager::
beginReadingLayer(const char *aLyrName,
                  pRPL::ReadingOption readOpt) {
  pRPL::Layer *pLyr = getLayerByName(aLyrName, true);
  if(pLyr == NULL) {
    _prc.abort();
    return false;
  }

  if(!pLyr->isDecomposed()) { // non-decomposed Layer
    if(readOpt == pRPL::PARA_READING) {
      if(!_prc.isMaster() && !_prc.isWriter()) {
        if(!pLyr->loadCellspaceByGDAL()) {
          return false;
        }
      }
    }
    else if(readOpt == pRPL::CENT_READING ||
            readOpt == pRPL::CENTDEL_READING) {
      if(_prc.isMaster()) {
        if(!pLyr->loadCellspaceByGDAL()) {
          _prc.abort();
          return false;
        }
      }
      if(!_iBcastCellspaceBegin(aLyrName)) {
        return false;
      }
    }
    else if(readOpt == pRPL::PGT_READING) {
      if(!_prc.isMaster() && !_prc.isWriter()) {
        if(!pLyr->loadCellspaceByPGTIOL()) {
          return false;
        }
      }
    }
  } // end -- non-decomposed Layer

  else { // decomposed Layer
    if(readOpt == pRPL::PARA_READING) { // parallel reading
      const pRPL::IntVect* pvSubspcIDs = _mOwnerships.subspcIDsOnPrc(_prc.id());
      if(pvSubspcIDs != NULL) {
        for(int iSub = 0; iSub < pvSubspcIDs->size(); iSub++) {
          int mappedSubspcID = pvSubspcIDs->at(iSub);
          if(!pLyr->loadSubCellspaceByGDAL(mappedSubspcID)) {
            return false;
          }
        }
      }
    } // end -- parallel reading

    else if(readOpt == pRPL::CENT_READING ||
            readOpt == pRPL::CENTDEL_READING) { // centralized reading
      if(_prc.isMaster()) { // master process
        pRPL::IntVect vMappedSubspcIDs = _mOwnerships.mappedSubspcIDs();
        for(int iSub = 0; iSub < vMappedSubspcIDs.size(); iSub++) {
          int mappedSubspcID = vMappedSubspcIDs[iSub];
          int prcID = _mOwnerships.findPrcBySubspcID(mappedSubspcID);
          if(mappedSubspcID >= 0 && prcID >= 0) {
            if(!pLyr->loadSubCellspaceByGDAL(mappedSubspcID)) {
              return false;
            }
            if(prcID != _prc.id()) { // transfer the SubCellspace if it belongs to a worker
              if(!_iTransfSubspaceBegin(_prc.id(), prcID, aLyrName, mappedSubspcID)) {
                return false;
              }
            }
          }
        } // end -- for(iSub)
      } // end -- master process

      else if(!_prc.isWriter()){ // worker process
        const pRPL::IntVect* pvSubspcIDs = _mOwnerships.subspcIDsOnPrc(_prc.id());
        if(pvSubspcIDs != NULL) {
          for(int iSub = 0; iSub < pvSubspcIDs->size(); iSub++) {
            int mappedSubspcID = pvSubspcIDs->at(iSub);
            if(!_iTransfSubspaceBegin(0, _prc.id(), aLyrName, mappedSubspcID)) {
              return false;
            }
          }
        }
      } // end -- worker process
    } // end -- centralized reading

    else if(readOpt == pRPL::PGT_READING) {
      const pRPL::IntVect* pvSubspcIDs = _mOwnerships.subspcIDsOnPrc(_prc.id());
      if(pvSubspcIDs != NULL) {
        for(int iSub = 0; iSub < pvSubspcIDs->size(); iSub++) {
          int mappedSubspcID = pvSubspcIDs->at(iSub);
          if(!pLyr->loadSubCellspaceByPGTIOL(mappedSubspcID)) {
            return false;
          }
        }
      }
    } // end -- pGTIOL reading

  } // end -- decomposed Layer
  return true;
}

void pRPL::DataManager::
finishReadingLayers(pRPL::ReadingOption readOpt) {
  // Complete all pending transfers, only for centralized I/O
  pRPL::IntVect vCmpltBcstIDs, vCmpltSendIDs, vCmpltRecvIDs;
  bool delSentSpc = (readOpt == pRPL::CENTDEL_READING) ? true : false;
  int nCmpltBcsts = 0, nCmpltSends = 0, nCmpltRecvs = 0;
  while(!_vBcstSpcReqs.empty() && nCmpltBcsts != MPI_UNDEFINED) {
    nCmpltBcsts = _iTransfSpaceTest(vCmpltBcstIDs, pRPL::BCST_TRANSF, delSentSpc, false, false);
  }
  while(!_vSendSpcReqs.empty() && nCmpltSends != MPI_UNDEFINED) {
    nCmpltSends = _iTransfSpaceTest(vCmpltSendIDs, pRPL::SEND_TRANSF, delSentSpc, false, false);
  }
  while(!_vRecvSpcReqs.empty() && nCmpltRecvs != MPI_UNDEFINED) {
    nCmpltRecvs = _iTransfSpaceTest(vCmpltRecvIDs, pRPL::RECV_TRANSF, delSentSpc, false, false);
  }

  _clearTransfSpaceReqs();
}

int pRPL::DataManager::
nNbrhds() const {
  return _lNbrhds.size();
}

pRPL::Neighborhood* pRPL::DataManager::
addNbrhd(const char *aNbrhdName) {
  pRPL::Neighborhood *pNbrhd = getNbrhdByName(aNbrhdName, false);
  if(pNbrhd != NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Warning: Neighborhood name (" << aNbrhdName \
         << ") already exists in DataManager. Failed to add new Neighborhood" \
         << endl;
  }
  else {
    _lNbrhds.push_back(pRPL::Neighborhood(aNbrhdName));
    pNbrhd = &(_lNbrhds.back());
  }
  return pNbrhd;
}

pRPL::Neighborhood* pRPL::DataManager::
addNbrhd(const pRPL::Neighborhood &nbrhd) {
  const char *aNbrhdName = nbrhd.name();
  pRPL::Neighborhood *pNbrhd = getNbrhdByName(aNbrhdName, false);
  if(pNbrhd != NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Warning: Neighborhood name (" << aNbrhdName \
         << ") already exists in DataManager. Failed to add new Neighborhood" \
        << endl;
  }
  else {
    _lNbrhds.push_back(pRPL::Neighborhood(nbrhd));
    pNbrhd = &(_lNbrhds.back());
  }
  return pNbrhd;
}

bool pRPL::DataManager::
rmvNbrhdByIdx(int nbrhdID,
              bool warning) {
  bool done = true;
  if(nbrhdID < 0 || nbrhdID >= _lNbrhds.size()) {
    done = false;
    if(warning) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: invalid Neighborhood index (" << nbrhdID \
           << ")" << endl;
    }
  }
  else {
    list<pRPL::Neighborhood>::iterator itrNbrhd = _lNbrhds.begin();
    advance(itrNbrhd, nbrhdID);
    _lNbrhds.erase(itrNbrhd);
  }
  return done;
}

bool pRPL::DataManager::
rmvNbrhdByName(const char *aNbrhdName,
               bool warning) {
  bool done = false;
  list<pRPL::Neighborhood>::iterator itrNbrhd = _lNbrhds.begin();
  while(itrNbrhd != _lNbrhds.end()) {
    if(strcmp(aNbrhdName, itrNbrhd->name()) == 0) {
      _lNbrhds.erase(itrNbrhd);
      done = true;
      break;
    }
    itrNbrhd++;
  }
  if(!done && warning) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: invalid Neighborhood name (" << aNbrhdName \
        << ")" << endl;
  }
  return done;
}

void pRPL::DataManager::
clearNbrhds() {
  _lNbrhds.clear();
}

pRPL::Neighborhood* pRPL::DataManager::
getNbrhdByIdx(int nbrhdID,
              bool warning) {
  pRPL::Neighborhood *pNbrhd = NULL;
  if(nbrhdID >= 0 && nbrhdID < _lNbrhds.size()) {
    list<pRPL::Neighborhood>::iterator itrNbrhd = _lNbrhds.begin();
    advance(itrNbrhd, nbrhdID);
    pNbrhd = &(*itrNbrhd);
  }
  else {
    if(warning) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Warning: NO Neighborhood with index (" << nbrhdID \
           << ") was found" << endl;
    }
  }
  return pNbrhd;
}

const pRPL::Neighborhood* pRPL::DataManager::
getNbrhdByIdx(int nbrhdID,
              bool warning) const {
  const pRPL::Neighborhood *pNbrhd = NULL;
  if(nbrhdID >= 0 && nbrhdID < _lNbrhds.size()) {
    list<pRPL::Neighborhood>::const_iterator itrNbrhd = _lNbrhds.begin();
    advance(itrNbrhd, nbrhdID);
    pNbrhd = &(*itrNbrhd);
  }
  else {
    if(warning) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Warning: NO Neighborhood with index (" << nbrhdID \
           << ") was found" << endl;
    }
  }
  return pNbrhd;
}

int pRPL::DataManager::
getNbrhdIdxByName(const char *aNbrhdName,
                  bool warning) {
  int nbrhdIdx = ERROR_ID, nbrhdCount = 0;
  bool found = false;
  list<pRPL::Neighborhood>::iterator itrNbrhd = _lNbrhds.begin();
  while(itrNbrhd != _lNbrhds.end()) {
    if(strcmp(aNbrhdName, itrNbrhd->name()) == 0) {
      nbrhdIdx = nbrhdCount;
      found = true;
      break;
    }
    nbrhdCount++;
    itrNbrhd++;
  }
  if(!found && warning) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Warning: NO Neighborhood with name (" \
         << aNbrhdName << ") was found" << endl;
  }
  return nbrhdIdx;
}

pRPL::Neighborhood* pRPL::DataManager::
getNbrhdByName(const char *aNbrhdName,
               bool warning) {
  pRPL::Neighborhood* pNbrhd = NULL;
  list<pRPL::Neighborhood>::iterator itrNbrhd = _lNbrhds.begin();
  while(itrNbrhd != _lNbrhds.end()) {
    if(strcmp(aNbrhdName, itrNbrhd->name()) == 0) {
      pNbrhd = &(*itrNbrhd);
      break;
    }
    itrNbrhd++;
  }
  if(pNbrhd == NULL && warning) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Warning: NO Neighborhood with name (" \
         << aNbrhdName << ") was found" << endl;
  }
  return pNbrhd;
}

const pRPL::Neighborhood* pRPL::DataManager::
getNbrhdByName(const char *aNbrhdName,
               bool warning) const {
  const pRPL::Neighborhood* pNbrhd = NULL;
  list<pRPL::Neighborhood>::const_iterator itrNbrhd = _lNbrhds.begin();
  while(itrNbrhd != _lNbrhds.end()) {
    if(strcmp(aNbrhdName, itrNbrhd->name()) == 0) {
      pNbrhd = &(*itrNbrhd);
      break;
    }
    itrNbrhd++;
  }
  if(pNbrhd == NULL && warning) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: NO Neighborhood with name (" \
        << aNbrhdName << ") was found" << endl;
  }
  return pNbrhd;
}

bool pRPL::DataManager::
dcmpLayer(int lyrID,
          int nbrhdID,
          int nRowSubspcs,
          int nColSubspcs) {
  pRPL::Layer *pLyr = getLayerByIdx(lyrID, true);
  if(pLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to apply decomposition on a NULL Layer" \
         << endl;
    return false;
  }

  pRPL::Neighborhood *pNbrhd = nbrhdID < 0 ? NULL : getNbrhdByIdx(nbrhdID, true);

  return pLyr->decompose(nRowSubspcs, nColSubspcs, pNbrhd);
}

bool pRPL::DataManager::
dcmpLayer(const char *aLyrName,
          const char *aNbrhdName,
          int nRowSubspcs,
          int nColSubspcs) {
  pRPL::Layer *pLyr = getLayerByName(aLyrName, true);
  if(pLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to apply decomposition on a NULL Layer" \
         << endl;
    return false;
  }

  pRPL::Neighborhood *pNbrhd = aNbrhdName == NULL ? NULL : getNbrhdByName(aNbrhdName, true);

  return pLyr->decompose(nRowSubspcs, nColSubspcs, pNbrhd);
}

bool pRPL::DataManager::
propagateDcmp(int fromLyrID,
              int toLyrID) {
  pRPL::Layer *pFromLyr = getLayerByIdx(fromLyrID, true);
  if(pFromLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to copy decomposition from a NULL Layer" \
         << endl;
    return false;
  }

  pRPL::Layer *pToLyr = getLayerByIdx(toLyrID, true);
  if(pToLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to copy decomposition to a NULL Layer" \
        << endl;
    return false;
  }

  if(pFromLyr == pToLyr) {
    return true;
  }

  return pToLyr->copyDcmp(pFromLyr);
}

bool pRPL::DataManager::
propagateDcmp(int fromLyrID,
              const IntVect &vToLyrIDs) {
  for(int iToLyrIdx = 0; iToLyrIdx < vToLyrIDs.size(); iToLyrIdx++) {
    int toLyrID = vToLyrIDs[iToLyrIdx];
    if(!propagateDcmp(fromLyrID, toLyrID)) {
      return false;
    }
  }
  return true;
}

bool pRPL::DataManager::
propagateDcmp(const string &fromLyrName,
              const string &toLyrName) {
  pRPL::Layer *pFromLyr = getLayerByName(fromLyrName.c_str(), true);
  if(pFromLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to copy decomposition from a NULL Layer" \
         << endl;
    return false;
  }

  pRPL::Layer *pToLyr = getLayerByName(toLyrName.c_str(), true);
  if(pToLyr == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to copy decomposition from a NULL Layer" \
        << endl;
    return false;
  }

  if(pFromLyr == pToLyr) {
    return true;
  }

  return pToLyr->copyDcmp(pFromLyr);
}

bool pRPL::DataManager::
propagateDcmp(const string &fromLyrName,
              const vector<string> &vToLyrNames) {
  for(int iToLyrName = 0; iToLyrName < vToLyrNames.size(); iToLyrName++) {
    const string &toLyrName = vToLyrNames[iToLyrName];
    if(!propagateDcmp(fromLyrName, toLyrName)) {
      return false;
    }
  }
  return true;
}

bool pRPL::DataManager::
dcmpLayers(const IntVect &vLyrIDs,
           int nbrhdID,
           int nRowSubspcs,
           int nColSubspcs) {
  if(!dcmpLayer(vLyrIDs[0], nbrhdID,
                nRowSubspcs, nColSubspcs)) {
    return false;
  }
  if(!propagateDcmp(vLyrIDs[0], vLyrIDs)) {
    return false;
  }
  return true;
}

bool pRPL::DataManager::
dcmpLayers(const vector<string> &vLyrNames,
           const char *aNbrhdName,
           int nRowSubspcs,
           int nColSubspcs) {
  if(!dcmpLayer(vLyrNames[0].c_str(), aNbrhdName,
                nRowSubspcs, nColSubspcs)) {
    return false;
  }
  if(!propagateDcmp(vLyrNames[0], vLyrNames)) {
    return false;
  }
  return true;
}

bool pRPL::DataManager::
dcmpLayers(pRPL::Transition &trans,
           int nRowSubspcs,
           int nColSubspcs) {
  const char *aNbrhdName = trans.getNbrhdName().empty() ? NULL : trans.getNbrhdName().c_str();
  if(!dcmpLayer(trans.getPrimeLyrName().c_str(),
                aNbrhdName,
                nRowSubspcs, nColSubspcs)) {
    return false;
  }

  if(!propagateDcmp(trans.getPrimeLyrName(), trans.getInLyrNames()) ||
     !propagateDcmp(trans.getPrimeLyrName(), trans.getOutLyrNames())) {
    return false;
  }
  return true;
}

bool pRPL::DataManager::
dcmpAllLayers(const char *aNbrhdName,
              int nRowSubspcs,
              int nColSubspcs) {
  list<pRPL::Layer>::iterator iLyr = _lLayers.begin();
  pRPL::Layer *p1stLyr = &(*iLyr);
  pRPL::Neighborhood *pNbrhd = aNbrhdName == NULL ? NULL : getNbrhdByName(aNbrhdName, true);
  if(!p1stLyr->decompose(nRowSubspcs, nColSubspcs, pNbrhd)) {
    return false;
  }
  iLyr++;

  while(iLyr != _lLayers.end()) {
    if(!(*iLyr).copyDcmp(p1stLyr)) {
      return false;
    }
    iLyr++;
  }

  return true;
}

const pRPL::OwnershipMap& pRPL::DataManager::
ownershipMap() const {
  return _mOwnerships;
}

bool pRPL::DataManager::
syncOwnershipMap() {
  vector<char> buf;
  int bufSize = 0;
  if(_prc.isMaster()) {
    if(!_mOwnerships.map2buf(buf)) {
      _prc.abort();
      return false;
    }
    bufSize = buf.size();
  }

  MPI_Bcast(&bufSize, 1, MPI_INT, 0, _prc.comm());

  if(!_prc.isMaster()) {
    buf.resize(bufSize);
  }

  MPI_Bcast(&(buf[0]), bufSize, MPI_CHAR, 0, _prc.comm());

  if(!_prc.isMaster()) {
    if(!_mOwnerships.buf2map(buf)) {
      _prc.abort();
      return false;
    }
  }

  return true;
}

bool pRPL::DataManager::
initTaskFarm(pRPL::Transition &trans,
             pRPL::MappingMethod mapMethod,
             int nSubspcs2Map,
             pRPL::ReadingOption readOpt) {
  if(_prc.nProcesses() < 2) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: unable to initialize task farming when there are " \
         << _prc.nProcesses() << " Processes" \
         << endl;
    _prc.abort();
    return false;
  }
  else if(_prc.nProcesses() <= 2 && _prc.hasWriter()) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to initialize task farming when there are " \
        << _prc.nProcesses() << " Processes that include a writer Process" \
        << endl;
    _prc.abort();
    return false;
  }

  pRPL::Layer *pPrimeLyr = getLayerByName(trans.getPrimeLyrName().c_str(), true);
  if(pPrimeLyr == NULL) {
    _prc.abort();
    return false;
  }
  if(!pPrimeLyr->isDecomposed()) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: the primary Layer (" << trans.getPrimeLyrName() \
         << ") is NOT decomposed for task farming" \
         << endl;
    _prc.abort();
    return false;
  }

  pRPL::IntVect vSubspcIDs = pPrimeLyr->allSubCellspaceIDs();
  pRPL::IntVect vPrcIDs = _prc.allPrcIDs(false, false); // exclude the master and writer processes
  _mOwnerships.init(vSubspcIDs, vPrcIDs);
  _mOwnerships.mapping(mapMethod, nSubspcs2Map); // initial mapping

  _clearTransfSpaceReqs();

  if(!_initReadInData(trans, readOpt)) {
    _prc.abort();
    return false;
  }

  return true;
}

pRPL::EvaluateReturn pRPL::DataManager::
evaluate_TF(pRPL::EvaluateType evalType,
            pRPL::Transition &trans,
            pRPL::ReadingOption readOpt,
            pRPL::WritingOption writeOpt,
            bool ifInitNoData,
            bool ifSyncOwnershipMap,
            const pRPL::LongVect *pvGlbIdxs) {
  if(trans.needExchange()) {
    if(!trans.onlyUpdtCtrCell() &&
        writeOpt != pRPL::NO_WRITING) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
          << " Error: when a non-centralized neighborhood-scope Transition is to be used," \
          << " and data exchange and writing are required, "
          << " only the static tasking can be used for load-balancing" \
          << endl;
      _prc.abort();
      return pRPL::EVAL_FAILED;
    }
    if(writeOpt == pRPL::CENTDEL_WRITING ||
       writeOpt == pRPL::PARADEL_WRITING) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: when data exchange is needed,"\
           << " deleting SubCellspaces during the writing is NOT allowed" \
           << endl;
      _prc.abort();
      return pRPL::EVAL_FAILED;
    }
  }

  if(trans.needExchange()) {
    if(!_calcExchangeBRs(trans)) {
      _prc.abort();
      return pRPL::EVAL_FAILED;
    }
  }

  if(_prc.isMaster() || _prc.isWriter()) {
    if(_prc.hasWriter()) {
      _nActiveWorkers = _prc.nProcesses() - 2; // exclude the master and writer processes
    }
    else {
      _nActiveWorkers = _prc.nProcesses() - 1; // exclude the master process
    }
  }

  pRPL::EvaluateReturn done;
  if(_prc.isMaster()) { // master process
    done = _masterTF(trans, readOpt, writeOpt);
  }
  else if(_prc.isWriter()) { // writer process
    done = _writerTF(trans, writeOpt);
  }
  else { // worker processes
    done = _workerTF(evalType, trans, readOpt, writeOpt, ifInitNoData, pvGlbIdxs);
  }

  if(done != pRPL::EVAL_FAILED) {
    if(trans.needExchange() || ifSyncOwnershipMap) {
      if(!syncOwnershipMap()) {
        done = pRPL::EVAL_FAILED;
      }
      else if(trans.needExchange()) {
        if(!_calcExchangeRoutes(trans)) {
          done = pRPL::EVAL_FAILED;
        }
        else if(!_makeCellStream(trans) ||
                !_iExchangeBegin() ||
                !_iExchangeEnd()) {
          done = pRPL::EVAL_FAILED;
        }
      }
    }
  }

  return done;
}

bool pRPL::DataManager::
initStaticTask(pRPL::Transition &trans,
               pRPL::MappingMethod mapMethod,
               pRPL::ReadingOption readOpt) {
  if(_prc.nProcesses() <= 1 && _prc.hasWriter()) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: unable to initialize static tasking when there are " \
        << _prc.nProcesses() << " Processes that include a writer Process" \
        << endl;
    _prc.abort();
    return false;
  }

  pRPL::Layer *pPrimeLyr = getLayerByName(trans.getPrimeLyrName().c_str(), true);
  if(pPrimeLyr == NULL) {
    _prc.abort();
    return false;
  }

  pRPL::IntVect vSubspcIDs = pPrimeLyr->allSubCellspaceIDs();
  pRPL::IntVect vPrcIDs = _prc.allPrcIDs(true, false); // exclude the writer processes
  _mOwnerships.init(vSubspcIDs, vPrcIDs);
  _mOwnerships.mapping(mapMethod); // map all SubCellspaces

  _clearTransfSpaceReqs();

  if(!_initReadInData(trans, readOpt)) {
    _prc.abort();
    return false;
  }

  // Complete all pending transfers, only for centralized I/O
  finishReadingLayers(readOpt);

  return true;
}

pRPL::EvaluateReturn pRPL::DataManager::
evaluate_ST(pRPL::EvaluateType evalType,
            pRPL::Transition &trans,
            pRPL::WritingOption writeOpt,
            bool ifInitNoData,
            const pRPL::LongVect *pvGlbIdxs) {
  if(!trans.onlyUpdtCtrCell() &&
     trans.needExchange() &&
     writeOpt != pRPL::NO_WRITING &&
     (_prc.hasWriter() ||
      (writeOpt != pRPL::CENT_WRITING || writeOpt != pRPL::CENTDEL_WRITING))) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: when a non-centralized neighborhood-scope Transition is to be used," \
        << " and data exchange and writing are required, " \
        << " only the centralized writing without a writer can be used" \
        << endl;
    _prc.abort();
    return pRPL::EVAL_FAILED;
  }

  if(trans.needExchange()) {
    if(!_calcExchangeBRs(trans) ||
       !_calcExchangeRoutes(trans)) {
      _prc.abort();
      return pRPL::EVAL_FAILED;
    }

    /*
      for(int prcID = 0; prcID < _prc.nProcesses(); prcID++) {
        if(_prc.id() == prcID) {
          cout << "Send Route on process " << prcID << endl;
          pRPL::ExchangeMap::iterator itrSendRoute = _mSendRoutes.begin();
          while(itrSendRoute != _mSendRoutes.end()) {
            cout << itrSendRoute->first << ": ";
            vector<pRPL::ExchangeNode> &vSendNodes = itrSendRoute->second;
            for(int iSend = 0; iSend < vSendNodes.size(); iSend++) {
              pRPL::ExchangeNode &sendNode = vSendNodes[iSend];
              cout << "(" << sendNode.subspcGlbID << ", " \
                  << sendNode.iDirection << ", "
                  << sendNode.iNeighbor << ")";
            }
            cout << endl;
            itrSendRoute++;
          }


          cout << "Recv Route on process " << prcID << endl;
          pRPL::ExchangeMap::iterator itrRecvRoute = _mRecvRoutes.begin();
          while(itrRecvRoute != _mRecvRoutes.end()) {
            cout << itrRecvRoute->first << ": ";
            vector<pRPL::ExchangeNode> &vRecvNodes = itrRecvRoute->second;
            for(int iRecv = 0; iRecv < vRecvNodes.size(); iRecv++) {
              pRPL::ExchangeNode &recvNode = vRecvNodes[iRecv];
              cout << "(" << recvNode.subspcGlbID << ", " \
                  << recvNode.iDirection << ", "
                  << recvNode.iNeighbor << ")";
            }
            cout << endl;
            itrRecvRoute++;
          }
        }
        _prc.sync();
      }
     */
  } // end -- if(trans.needExchange())
  if(_prc.isMaster()) {
    if(_prc.hasWriter()) {
      _nActiveWorkers = _prc.nProcesses() - 2; // exclude the master and writer processes
    }
    else {
      _nActiveWorkers = _prc.nProcesses() - 1; // exclude the master process
    }
  }
  else if(_prc.isWriter()) {
    _nActiveWorkers = _prc.nProcesses() - 1; // exclude the writer process
  }
  pRPL::EvaluateReturn done;
  if(_prc.isWriter()) { // writer process
    done = _writerST(trans, writeOpt);
  }
  else { // worker processes
    done = _workerST(evalType, trans, writeOpt, ifInitNoData, pvGlbIdxs);
    if(_prc.isMaster() && done != pRPL::EVAL_FAILED) { // master process
      done = _masterST(trans, writeOpt);
    }
  }
  return done;
}

bool pRPL::DataManager::
mergeTmpSubspcByGDAL(pRPL::Transition &trans,
                     int subspcGlbID) {
  if((_prc.hasWriter() && _prc.isWriter()) ||
     (!_prc.hasWriter() && _prc.isMaster())) {
    const vector<string> &vOutLyrs = trans.getOutLyrNames();
    for(int iOutLyr = 0; iOutLyr < vOutLyrs.size(); iOutLyr++) {
      pRPL::Layer *pOutLyr = getLayerByName(vOutLyrs[iOutLyr].c_str(), true);
      if(pOutLyr == NULL) {
        _prc.abort();
        return false;
      }
      if(!pOutLyr->mergeSubCellspaceByGDAL(subspcGlbID, _prc.nProcesses())) {
        _prc.abort();
        return false;
      }
      if(!_mOwnerships.findSubspcIDOnPrc(subspcGlbID, _prc.id())) {
        pOutLyr->delSubCellspace_glbID(subspcGlbID);
      }
    }
  }
  return true;
}

bool pRPL::DataManager::
mergeAllTmpSubspcsByGDAL(pRPL::Transition &trans) {
  if((_prc.hasWriter() && _prc.isWriter()) ||
     (!_prc.hasWriter() && _prc.isMaster())) {
    const vector<string> &vOutLyrs = trans.getOutLyrNames();
    for(int iOutLyr = 0; iOutLyr < vOutLyrs.size(); iOutLyr++) {
      pRPL::Layer *pOutLyr = getLayerByName(vOutLyrs[iOutLyr].c_str(), true);
      if(pOutLyr == NULL) {
        _prc.abort();
        return false;
      }
      list<pRPL::SubCellspaceInfo>::const_iterator itrSubInfo = pOutLyr->allSubCellspaceInfos()->begin();
      while(itrSubInfo != pOutLyr->allSubCellspaceInfos()->end()) {
        int glbID = itrSubInfo->id();
        if(!pOutLyr->mergeSubCellspaceByGDAL(glbID, _prc.nProcesses())) {
          return false;
        }
        if(!_mOwnerships.findSubspcIDOnPrc(glbID, _prc.id())) {
          pOutLyr->delSubCellspace_glbID(glbID);
        }
        itrSubInfo++;
      }
    } // end of iOutLyr loop
  }
  return true;
}
