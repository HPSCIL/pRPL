#include "prpl-transition.h"

pRPL::Transition::
Transition(bool onlyUpdtCtrCell,
           bool needExchange,
           bool edgesFirst)
    :_pNbrhd(NULL),
     _onlyUpdtCtrCell(onlyUpdtCtrCell),
     _needExchange(needExchange),
     _edgesFirst(edgesFirst) {}

void pRPL::Transition::
addInputLyr(const char *aInLyrName,
            bool isPrimeLyr) {
  if(aInLyrName != NULL &&
     std::find(_vInLyrNames.begin(), _vInLyrNames.end(), aInLyrName) == _vInLyrNames.end()) {
    _vInLyrNames.push_back(aInLyrName);
    if(_mpCellspcs.find(aInLyrName) == _mpCellspcs.end()) {
      _mpCellspcs[aInLyrName] = NULL;
    }
    if(isPrimeLyr) {
      _primeLyrName = aInLyrName;
    }
  }
}

void pRPL::Transition::
addOutputLyr(const char *aOutLyrName,
             bool isPrimeLyr) {
  if(aOutLyrName != NULL &&
     std::find(_vOutLyrNames.begin(), _vOutLyrNames.end(), aOutLyrName) == _vOutLyrNames.end()) {
    _vOutLyrNames.push_back(aOutLyrName);
    if(_mpCellspcs.find(aOutLyrName) == _mpCellspcs.end()) {
      _mpCellspcs[aOutLyrName] = NULL;
    }
    if(isPrimeLyr) {
      _primeLyrName = aOutLyrName;
    }
  }
}

bool pRPL::Transition::
setLyrsByNames(vector<string> *pvInLyrNames,
               vector<string> *pvOutLyrNames,
               string &primeLyrName) {
  clearLyrSettings();

  vector<string>::iterator itrLyrName;
  bool foundPrm = false;
  if(pvInLyrNames != NULL) {
    itrLyrName = std::find(pvInLyrNames->begin(), pvInLyrNames->end(), primeLyrName);
    if(itrLyrName != pvInLyrNames->end()) {
      _primeLyrName = primeLyrName;
      foundPrm = true;
    }
    _vInLyrNames = *pvInLyrNames;
  }
  if(pvOutLyrNames != NULL) {
    itrLyrName = find(pvOutLyrNames->begin(), pvOutLyrNames->end(), primeLyrName);
    if(itrLyrName != pvOutLyrNames->end()) {
      _primeLyrName = primeLyrName;
      foundPrm = true;
    }
    _vOutLyrNames = *pvOutLyrNames;
  }
  if(!foundPrm) {
    cerr << __FILE__ << " function: " << __FUNCTION__ \
         << " Error: primary Layer name (" << primeLyrName \
         << ") is not found in either Input Layer names or Output Layer names" << endl;
    return false;
  }

  for(int iLyrName = 0; iLyrName < _vInLyrNames.size(); iLyrName++) {
    const string &lyrName = _vInLyrNames[iLyrName];
    _mpCellspcs[lyrName] = NULL;
  }
  for(int iLyrName = 0; iLyrName < _vOutLyrNames.size(); iLyrName++) {
    const string &lyrName = _vOutLyrNames[iLyrName];
    if(std::find(_vInLyrNames.begin(), _vInLyrNames.end(), lyrName) != _vInLyrNames.end()) {
      /*
      cerr << __FILE__ << " function: " << __FUNCTION__ \
           << " Error: output Layer (" << lyrName \
           << ") can NOT be also input Layer" << endl;
      return false;
      */
      continue; // ignore the output Layer that is also an input Layer
    }
    _mpCellspcs[lyrName] = NULL;
  }

  return true;
}

const vector<string>& pRPL::Transition::
getInLyrNames() const {
  return _vInLyrNames;
}

const vector<string>& pRPL::Transition::
getOutLyrNames() const {
  return _vOutLyrNames;
}

const string& pRPL::Transition::
getPrimeLyrName() const {
  return _primeLyrName;
}

bool pRPL::Transition::
isInLyr(const string &lyrName) const {
  return (std::find(_vInLyrNames.begin(), _vInLyrNames.end(), lyrName) != _vInLyrNames.end());
}

bool pRPL::Transition::
isOutLyr(const string &lyrName) const {
  return (std::find(_vOutLyrNames.begin(), _vOutLyrNames.end(), lyrName) != _vOutLyrNames.end());
}

bool pRPL::Transition::
isPrimeLyr(const string &lyrName) const {
  return (lyrName == _primeLyrName);
}

void pRPL::Transition::
clearLyrSettings() {
  _vInLyrNames.clear();
  _vOutLyrNames.clear();
  _primeLyrName.clear();
  _mpCellspcs.clear();
}

bool pRPL::Transition::
setCellspace(const string &lyrName,
             pRPL::Cellspace *pCellspc) {
  map<string, pRPL::Cellspace *>::iterator itrCellspc = _mpCellspcs.find(lyrName);
  if(itrCellspc == _mpCellspcs.end()) {
    cerr << __FILE__ << " function: " << __FUNCTION__ \
         << " Error: unable to find a Layer with the name (" \
         << lyrName << ")" << endl;
    return false;
  }

  if(pCellspc == NULL ||
     pCellspc->isEmpty(true)) {
    cerr << __FILE__ << " function: " << __FUNCTION__ \
         << " Error: NULL or empty Cellspace to be added to the Transition" \
         << endl;
    return false;
  }
  itrCellspc->second = pCellspc;

  return true;
}

void pRPL::Transition::
clearCellspaces() {
  map<string, pRPL::Cellspace *>::iterator itrCellspc = _mpCellspcs.begin();
  while(itrCellspc != _mpCellspcs.end()) {
    itrCellspc->second = NULL;
    itrCellspc++;
  }
}

pRPL::Cellspace* pRPL::Transition::
getCellspaceByLyrName(const string &lyrName) {
  pRPL::Cellspace *pCellspc = NULL;
  map<string, pRPL::Cellspace *>::iterator itrCellspc = _mpCellspcs.find(lyrName);
  if(itrCellspc != _mpCellspcs.end()) {
    pCellspc = itrCellspc->second;
  }
  return pCellspc;
}

const pRPL::Cellspace* pRPL::Transition::
getCellspaceByLyrName(const string &lyrName) const {
  const pRPL::Cellspace *pCellspc = NULL;
  map<string, pRPL::Cellspace *>::const_iterator itrCellspc = _mpCellspcs.find(lyrName);
  if(itrCellspc != _mpCellspcs.end()) {
    pCellspc = itrCellspc->second;
  }
  return pCellspc;
}

void pRPL::Transition::
setUpdateTracking(bool toTrack) {
  for(int iOutLyrName = 0; iOutLyrName < _vOutLyrNames.size(); iOutLyrName++) {
    map<string, Cellspace *>::iterator itrCellspc = _mpCellspcs.find(_vOutLyrNames[iOutLyrName]);
    if(itrCellspc != _mpCellspcs.end() &&
       itrCellspc->second != NULL) {
      itrCellspc->second->setUpdateTracking(toTrack);
    }
  }
}

void pRPL::Transition::
clearUpdateTracks() {
  for(int iOutLyrName = 0; iOutLyrName < _vOutLyrNames.size(); iOutLyrName++) {
    map<string, Cellspace *>::iterator itrCellspc = _mpCellspcs.find(_vOutLyrNames[iOutLyrName]);
    if(itrCellspc != _mpCellspcs.end() &&
       itrCellspc->second != NULL) {
      itrCellspc->second->clearUpdatedIdxs();
    }
  }
}

void pRPL::Transition::
setNbrhdByName(const char *aNbrhdName) {
  if(aNbrhdName != NULL) {
    _nbrhdName = aNbrhdName;
  }
}

const string& pRPL::Transition::
getNbrhdName() const {
  return _nbrhdName;
}

void pRPL::Transition::
clearNbrhdSetting() {
  _nbrhdName.clear();
  _pNbrhd = NULL;
}

void pRPL::Transition::
setNbrhd(Neighborhood *pNbrhd) {
  _pNbrhd = pNbrhd;
}

pRPL::Neighborhood* pRPL::Transition::
getNbrhd() {
  return _pNbrhd;
}

const pRPL::Neighborhood* pRPL::Transition::
getNbrhd() const {
  return _pNbrhd;
}

void pRPL::Transition::
clearDataSettings() {
  clearLyrSettings();
  clearNbrhdSetting();
}

bool pRPL::Transition::
onlyUpdtCtrCell() const {
  return _onlyUpdtCtrCell;
}

bool pRPL::Transition::
needExchange() const {
  return _needExchange;
}

bool pRPL::Transition::
edgesFirst() const {
  return _edgesFirst;
}

pRPL::EvaluateReturn pRPL::Transition::
evalBR(const pRPL::CoordBR &br) {
  pRPL::Cellspace *pPrmSpc = getCellspaceByLyrName(getPrimeLyrName());
  if(pPrmSpc == NULL) {
    cerr << __FILE__ << " function: " << __FUNCTION__ \
         << " Error: unable to find the primary Cellspace with name (" \
         << getPrimeLyrName() << ")" << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::EvaluateReturn done = pRPL::EVAL_SUCCEEDED;
  if(br.valid(pPrmSpc->info()->dims())) {
    for(long iRow = br.minIRow(); iRow <= br.maxIRow(); iRow++) {
      for(long iCol = br.minICol(); iCol <= br.maxICol(); iCol++) {
        done = evaluate(pRPL::CellCoord(iRow, iCol));
        if(done == pRPL::EVAL_FAILED ||
           done == pRPL::EVAL_TERMINATED) {
          return done;
        }
      }
    }
  }

  return pRPL::EVAL_SUCCEEDED;
}

pRPL::EvaluateReturn pRPL::Transition::
evalRandomly(const pRPL::CoordBR &br) {
  pRPL::Cellspace *pPrmSpc = getCellspaceByLyrName(getPrimeLyrName());
  if(pPrmSpc == NULL) {
    cerr << __FILE__ << " function: " << __FUNCTION__ \
        << " Error: unable to find the primary Cellspace with name (" \
        << getPrimeLyrName() << ")" << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::EvaluateReturn done = pRPL::EVAL_SUCCEEDED;
  if(br.valid(pPrmSpc->info()->dims())) {
    while(done != pRPL::EVAL_TERMINATED && done != pRPL::EVAL_FAILED) {
      long iRow = rand() % br.nRows() + br.minIRow();
      long iCol = rand() % br.nCols() + br.minICol();
      pRPL::CellCoord coord2Eval(iRow, iCol);
      if(br.ifContain(coord2Eval)) {
        done = evaluate(coord2Eval);
      }
    }
  }

  return done;
}

pRPL::EvaluateReturn pRPL::Transition::
evalSelected(const pRPL::CoordBR &br,
             const pRPL::LongVect &vlclIdxs) {
  pRPL::Cellspace *pPrmSpc = getCellspaceByLyrName(getPrimeLyrName());
  if(pPrmSpc == NULL) {
    cerr << __FILE__ << " function: " << __FUNCTION__ \
        << " Error: unable to find the primary Cellspace with name (" \
        << getPrimeLyrName() << ")" << endl;
    return pRPL::EVAL_FAILED;
  }

  pRPL::EvaluateReturn done = pRPL::EVAL_SUCCEEDED;
  if(br.valid(pPrmSpc->info()->dims())) {
    for(int iIdx = 0; iIdx < vlclIdxs.size(); iIdx++) {
      pRPL::CellCoord coord = pPrmSpc->info()->idx2coord(vlclIdxs[iIdx]);
      if(br.ifContain(coord)) {
        done = evaluate(coord);
        if(done == pRPL::EVAL_FAILED ||
           done == pRPL::EVAL_TERMINATED) {
          return done;
        }
      }
    }
  }

  return pRPL::EVAL_SUCCEEDED;
}

bool pRPL::Transition::
afterSetCellspaces(int subCellspcGlbIdx) {
  return true;
}

bool pRPL::Transition::
afterSetNbrhd() {
  return true;
}

bool pRPL::Transition::
check() const {
  return true;
}

pRPL::EvaluateReturn pRPL::Transition::
evaluate(const CellCoord &coord) {
  return pRPL::EVAL_SUCCEEDED;
}

double pRPL::Transition::
workload(const CoordBR &workBR) const {
  return 0.0;
}
