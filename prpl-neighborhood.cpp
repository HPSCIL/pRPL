#include "prpl-neighborhood.h"

//#include <sys/_types/_size_t.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

/****************************************
*             Private Methods           *
*****************************************/
bool pRPL::Neighborhood::
_checkNbrCoords(const vector<pRPL::CellCoord> &vNbrCoords) const {
  bool valid = true;
  if(vNbrCoords.empty()) {
    cerr << "Error: the coordinate vector is empty" << endl;
    valid = false;
  }
  else if(std::find(vNbrCoords.begin(), vNbrCoords.end(), CellCoord(0, 0)) == vNbrCoords.end()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the coordinate vector has to have a coordinate [0, 0]" << endl;
    valid = false;
  }
  return valid;
}

bool pRPL::Neighborhood::
_validNumDirects() const {
  bool numDirectsValid = true;
  if(_mNbrIDMap.size() != 8) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid _mNbrIDMap size" << endl;
    numDirectsValid = false;
  }
  return numDirectsValid;
}

/****************************************
*             Public Methods            *
*****************************************/
pRPL::Neighborhood::
Neighborhood(const char *aName,
             pRPL::EdgeOption edgeOption,
             double virtualEdgeVal) {
  if(aName != NULL) {
    _name.assign(aName);
  }
  else {
    _name.assign("Untitled_Neighborhood");
  }
  _edgeOption = edgeOption;
  _pVirtualEdgeVal = NULL;

  if(edgeOption == pRPL::CUSTOM_VIRTUAL_EDGES) {
    _pVirtualEdgeVal = new double;
    memcpy(_pVirtualEdgeVal, &virtualEdgeVal, sizeof(double));
  }
}

pRPL::Neighborhood::
Neighborhood(const pRPL::Neighborhood &rhs) {
  _name = rhs._name;
  _vNbrs = rhs._vNbrs;
  _MBR = rhs._MBR;
  _mNbrIDMap = rhs._mNbrIDMap;
  _edgeOption = rhs._edgeOption;
  if(rhs._pVirtualEdgeVal != NULL) {
    _pVirtualEdgeVal = new double(*(rhs._pVirtualEdgeVal));
  }
  else {
    _pVirtualEdgeVal = NULL;
  }
}

pRPL::Neighborhood::
~Neighborhood() {
  if(_pVirtualEdgeVal != NULL) {
    delete _pVirtualEdgeVal;
  }
}

bool pRPL::Neighborhood::
init(const vector<pRPL::CellCoord> &vNbrCoords,
     double weight,
     pRPL::EdgeOption edgeOption,
     double virtualEdgeVal) {
  if(!_checkNbrCoords(vNbrCoords)) {
    return false;
  }

  _vNbrs.clear();
  _mNbrIDMap.clear();
  _mNbrIDMap.resize(8);

  _edgeOption = edgeOption;
  _pVirtualEdgeVal = NULL;
  if(edgeOption == pRPL::CUSTOM_VIRTUAL_EDGES) {
    _pVirtualEdgeVal = new double;
    memcpy(_pVirtualEdgeVal, &virtualEdgeVal, sizeof(double));
  }

  CellCoord nwCorner, seCorner;
  for(int iNbr = 0; iNbr < vNbrCoords.size(); iNbr++) {
    vector<pRPL::WeightedCellCoord>::iterator nbrItr = _vNbrs.begin();
    while(nbrItr != _vNbrs.end()) {
      pRPL::CellCoord coord(nbrItr->iRow(), nbrItr->iCol());
      if(coord == vNbrCoords[iNbr]) {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: a duplicated coordinate [" \
             << vNbrCoords[iNbr] \
             << "] has been found" << endl;
        return false;
      }
      nbrItr++;
    }

    _vNbrs.push_back(pRPL::WeightedCellCoord(vNbrCoords[iNbr].iRow(), vNbrCoords[iNbr].iCol(), weight));

    // calc _MBR
    if(0 == iNbr) {
      nwCorner = vNbrCoords[iNbr];
      seCorner = vNbrCoords[iNbr];
    }
    else{
      if(nwCorner.iRow() > vNbrCoords[iNbr].iRow()) {
         nwCorner.iRow(vNbrCoords[iNbr].iRow());
      }
      if(nwCorner.iCol() > vNbrCoords[iNbr].iCol()) {
         nwCorner.iCol(vNbrCoords[iNbr].iCol());
      }
      if(seCorner.iRow() < vNbrCoords[iNbr].iRow()) {
         seCorner.iRow(vNbrCoords[iNbr].iRow());
      }
      if(seCorner.iCol() < vNbrCoords[iNbr].iCol()) {
         seCorner.iCol(vNbrCoords[iNbr].iCol());
      }
    }

    // calc _mNbrIDMap
    if(vNbrCoords[iNbr].iRow() < 0) {
      if(vNbrCoords[iNbr].iCol() < 0) {
        _mNbrIDMap[pRPL::NORTHWEST_DIR].push_back(iNbr);
      }
      else if(vNbrCoords[iNbr].iCol() > 0) {
        _mNbrIDMap[pRPL::NORTHEAST_DIR].push_back(iNbr);
      }
      _mNbrIDMap[pRPL::NORTH_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iRow() > 0) {
      if(vNbrCoords[iNbr].iCol() < 0) {
        _mNbrIDMap[pRPL::SOUTHWEST_DIR].push_back(iNbr);
      }
      else if(vNbrCoords[iNbr].iCol() > 0) {
        _mNbrIDMap[pRPL::SOUTHEAST_DIR].push_back(iNbr);
      }
      _mNbrIDMap[pRPL::SOUTH_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iCol() < 0) {
      _mNbrIDMap[pRPL::WEST_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iCol() > 0) {
      _mNbrIDMap[pRPL::EAST_DIR].push_back(iNbr);
    }
  } // End of iNbr loop

  _MBR.nwCorner(nwCorner);
  _MBR.seCorner(seCorner);

  return true;
}

bool pRPL::Neighborhood::
init(const vector<pRPL::WeightedCellCoord> &vNbrCoords,
     pRPL::EdgeOption edgeOption,
     double virtualEdgeVal) {
  if(!_checkNbrCoords((const vector<pRPL::CellCoord> &)vNbrCoords)) {
    return false;
  }
  
  _vNbrs.clear();
  _mNbrIDMap.clear();
  _mNbrIDMap.resize(8);
  
  _edgeOption = edgeOption;
  _pVirtualEdgeVal = NULL;
  if(edgeOption == pRPL::CUSTOM_VIRTUAL_EDGES) {
    _pVirtualEdgeVal = new double;
    memcpy(_pVirtualEdgeVal, &virtualEdgeVal, sizeof(double));
  }

  pRPL::CellCoord nwCorner, seCorner;
  for(int iNbr = 0; iNbr < vNbrCoords.size(); iNbr++) {
    if(std::find(_vNbrs.begin(), _vNbrs.end(), vNbrCoords[iNbr]) != _vNbrs.end()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
      << " Error: a duplicated coordinate [" \
      << vNbrCoords[iNbr] \
      << "] has been found" << endl;
      return false;
    }
    
    _vNbrs.push_back(pRPL::WeightedCellCoord(vNbrCoords[iNbr]));
    
    // calc _MBR
    if(0 == iNbr) {
      nwCorner = vNbrCoords[iNbr];
      seCorner = vNbrCoords[iNbr];
    }
    else{
      if(nwCorner.iRow() > vNbrCoords[iNbr].iRow()) {
        nwCorner.iRow(vNbrCoords[iNbr].iRow());
      }
      if(nwCorner.iCol() > vNbrCoords[iNbr].iCol()) {
        nwCorner.iCol(vNbrCoords[iNbr].iCol());
      }
      if(seCorner.iRow() < vNbrCoords[iNbr].iRow()) {
        seCorner.iRow(vNbrCoords[iNbr].iRow());
      }
      if(seCorner.iCol() < vNbrCoords[iNbr].iCol()) {
        seCorner.iCol(vNbrCoords[iNbr].iCol());
      }
    }
    
    // calc _mNbrIDMap
    if(vNbrCoords[iNbr].iRow() < 0) {
      if(vNbrCoords[iNbr].iCol() < 0) {
        _mNbrIDMap[pRPL::NORTHWEST_DIR].push_back(iNbr);
      }
      else if(vNbrCoords[iNbr].iCol() > 0) {
        _mNbrIDMap[pRPL::NORTHEAST_DIR].push_back(iNbr);
      }
      _mNbrIDMap[pRPL::NORTH_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iRow() > 0) {
      if(vNbrCoords[iNbr].iCol() < 0) {
        _mNbrIDMap[pRPL::SOUTHWEST_DIR].push_back(iNbr);
      }
      else if(vNbrCoords[iNbr].iCol() > 0) {
        _mNbrIDMap[pRPL::SOUTHEAST_DIR].push_back(iNbr);
      }
      _mNbrIDMap[pRPL::SOUTH_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iCol() < 0) {
      _mNbrIDMap[pRPL::WEST_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iCol() > 0) {
      _mNbrIDMap[pRPL::EAST_DIR].push_back(iNbr);
    }
  } // End of iNbr loop
  
  _MBR.nwCorner(nwCorner);
  _MBR.seCorner(seCorner);
  
  return true;
}

bool pRPL::Neighborhood::
initSingleCell(double weight) {
  vector<pRPL::CellCoord> vNbrCoords;
  vNbrCoords.push_back(pRPL::CellCoord(0, 0));
  return init(vNbrCoords, weight, pRPL::FORBID_VIRTUAL_EDGES);
}

bool pRPL::Neighborhood::
initVonNeumann(double weight,
               pRPL::EdgeOption edgeOption,
               double virtualEdgeVal) {
  vector<pRPL::CellCoord> vNbrCoords;
  vNbrCoords.push_back(pRPL::CellCoord(0, 0));
  vNbrCoords.push_back(pRPL::CellCoord(-1, 0));
  vNbrCoords.push_back(pRPL::CellCoord(0, 1));
  vNbrCoords.push_back(pRPL::CellCoord(1, 0));
  vNbrCoords.push_back(pRPL::CellCoord(0, -1));
  return init(vNbrCoords, weight, edgeOption, virtualEdgeVal);
}

bool pRPL::Neighborhood::
initMoore(long nEdgeLength,
          double weight,
          pRPL::EdgeOption edgeOption,
          double virtualEdgeVal) {
  if(nEdgeLength <= 0 || nEdgeLength % 2 == 0) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid edge length (" << nEdgeLength \
         << "). A Moore neighborhood's edge length must be a positive odd number" << endl;
    return false;
  }

  vector<pRPL::CellCoord> vNbrCoords;
  vNbrCoords.push_back(pRPL::CellCoord(0, 0));
  for(long iRow = -(nEdgeLength/2); iRow <= nEdgeLength/2; iRow++) {
    for(long iCol = -(nEdgeLength/2); iCol <= nEdgeLength/2; iCol++) {
      if(iRow == 0 && iCol == 0) {
        continue;
      }
      vNbrCoords.push_back(pRPL::CellCoord(iRow, iCol));
    }
  }

  return init(vNbrCoords, weight, edgeOption, virtualEdgeVal);
}

void pRPL::Neighborhood::
clear() {
  _name.clear();
  _vNbrs.clear();
  _mNbrIDMap.clear();
  _MBR = CoordBR();
  _edgeOption = pRPL::FORBID_VIRTUAL_EDGES;
  if(_pVirtualEdgeVal != NULL) {
    delete _pVirtualEdgeVal;
    _pVirtualEdgeVal = NULL;
  }
}

pRPL::WeightedCellCoord& pRPL::Neighborhood::
operator[](int iNbr) {
  return _vNbrs.at(iNbr);
}

const pRPL::WeightedCellCoord& pRPL::Neighborhood::
operator[](int iNbr) const {
  return _vNbrs.at(iNbr);
}

pRPL::Neighborhood& pRPL::Neighborhood::
operator=(const pRPL::Neighborhood &rhs) {
  if(this != &rhs) {
    clear();
    _vNbrs = rhs._vNbrs;
    _MBR = rhs._MBR;
    _mNbrIDMap = rhs._mNbrIDMap;
    _edgeOption = rhs._edgeOption;
    if(rhs._pVirtualEdgeVal != NULL) {
      _pVirtualEdgeVal = new double(*(rhs._pVirtualEdgeVal));
    }
  }
  return *this;
}

const char* pRPL::Neighborhood::
name() const {
  return _name.c_str();
}

void pRPL::Neighborhood::
name(const char *aName) {
  if(aName != NULL) {
    _name.assign(aName);
  }
}

bool pRPL::Neighborhood::
isEmpty() const{
  return _vNbrs.empty();
}

int pRPL::Neighborhood::
size() const{
  return _vNbrs.size();  
}

bool pRPL::Neighborhood::
isEquallyWeighted(double &weight) const {
  bool equal = true;
  weight = 0.0;
  if(_vNbrs.empty()) {
    return equal;
  }
  weight = _vNbrs[0].weight();
  for(int iNbr = 1; iNbr < _vNbrs.size(); iNbr++) {
    if(weight != _vNbrs[iNbr].weight()) {
      equal = false;
      break;
    }
  }
  return equal;
}

long pRPL::Neighborhood::
minIRow() const {
  return _MBR.minIRow();
}

long pRPL::Neighborhood::
minICol() const {
  return _MBR.minICol();
}

long pRPL::Neighborhood::
maxIRow() const {
  return _MBR.maxIRow();
}

long pRPL::Neighborhood::
maxICol() const {
  return _MBR.maxICol();
}

long pRPL::Neighborhood::
nRows() const {
  return _MBR.nRows();
}

long pRPL::Neighborhood::
nCols() const {
  return _MBR.nCols();
}

const pRPL::CoordBR* pRPL::Neighborhood::
getMBR() const {
  return &(_MBR);
}

pRPL::EdgeOption pRPL::Neighborhood::
edgeOption() const {
  return _edgeOption;
}

const double* pRPL::Neighborhood::
virtualEdgeVal() const {
  return _pVirtualEdgeVal;
}

bool pRPL::Neighborhood::
virtualEdgeVal(const double veVal) {
  if(_edgeOption != pRPL::CUSTOM_VIRTUAL_EDGES) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: can NOT customize the virtual edge value," \
         << " because the Neighborhood\'s edge option is NOT CUSTOM_VIRTUAL_EDGE" \
         << endl;
    return false;
  }
  if(_pVirtualEdgeVal == NULL) {
    _pVirtualEdgeVal = new double(veVal);
  }
  else {
    memcpy(_pVirtualEdgeVal, &veVal, sizeof(double));
  }
  return true;
}

bool pRPL::Neighborhood::
hasNbrs(pRPL::MeshDir dir) const {
  bool nbrExist;
  if(!_validNumDirects()) {
    nbrExist = false; 
  }
  else {
    if(!_mNbrIDMap[dir].empty()) {
      nbrExist = true;
    }
    else {
      nbrExist = false;
    }
  }
  return nbrExist;
}

const pRPL::IntVect* pRPL::Neighborhood::
nbrIDs(pRPL::MeshDir dir) const {
  return &(_mNbrIDMap.at(dir));
}

bool pRPL::Neighborhood::
calcWorkBR(pRPL::CoordBR &workBR,
           const pRPL::SpaceDims &dims) const {
  if(isEmpty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to calculate workBR" \
         << " using an empty Neighborhood" << endl;
    return false;
  }
  if(dims.isNone()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to calculate workBR" \
         << " using an empty SpaceDims" << endl;
    return false;
  }
  if(dims.nRows() < nRows() ||
     dims.nCols() < nCols()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the dimensions (" << dims << ") is smaller than" \
         << " the dimensions of the Neighborhood (" << nRows() \
         << " " << nCols() << ")" << endl;
    return false;
  }

  if(_edgeOption == FORBID_VIRTUAL_EDGES) {
    workBR.nwCorner(-minIRow(), -minICol());
    workBR.seCorner(dims.nRows() - maxIRow() - 1,
                    dims.nCols() - maxICol() - 1);

    /*
    if((maxIRow() > 0 && workBR.nRows() < maxIRow()) ||
       (minIRow() < 0 && workBR.nRows() < -minIRow()) ||
       (maxICol() > 0 && workBR.nCols() < maxICol()) ||
       (minICol() < 0 && workBR.nCols() < -minICol())) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: workBR (" << workBR \
           << ") is too small to accommodate the Neighborhood" << endl;
      return false;
    }
    */
  }
  else {
    workBR.nwCorner(0, 0);
    workBR.seCorner(dims.nRows() - 1,
                    dims.nCols() - 1);
  }

  if(!workBR.valid(dims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid workBR (" << workBR \
         << ") has been produced" << endl;
    return false;
  }

  return true;
}

void pRPL::Neighborhood::
add2Buf(pRPL::CharVect &vBuf) const {
  size_t optionSize = sizeof(pRPL::EdgeOption);
  size_t coordSize = sizeof(pRPL::WeightedCellCoord);
  size_t intSize = sizeof(int);

  int bufSize = vBuf.size();
  vBuf.resize(bufSize + optionSize + intSize + size()*coordSize);

  memcpy(&(vBuf[bufSize]), &(_edgeOption), optionSize);
  bufSize += optionSize;

  int nNbrs = size();
  memcpy(&(vBuf[bufSize]), &(nNbrs), intSize);
  bufSize += intSize;

  for(int iNbr = 0; iNbr < size(); iNbr++) {
    memcpy(&(vBuf[bufSize]), &(_vNbrs[iNbr]), coordSize);
    bufSize += coordSize;
  }
}

bool pRPL::Neighborhood::
fromBuf(const pRPL::CharVect &vBuf,
        int &iChar) {
  size_t optionSize = sizeof(pRPL::EdgeOption);
  size_t coordSize = sizeof(pRPL::WeightedCellCoord);
  size_t intSize = sizeof(int);

  if(vBuf.size() - iChar < optionSize) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: buffer is too small to extract a Neighborhood" << endl;
    return false;
  }
  pRPL::EdgeOption edgeOption;
  memcpy(&edgeOption, &(vBuf[iChar]), optionSize);
  iChar += optionSize;

  if(vBuf.size() - iChar < intSize) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: buffer is too small to extract a Neighborhood" << endl;
    return false;
  }
  int nNbrs;
  memcpy(&nNbrs, &(vBuf[iChar]), intSize);
  iChar += intSize;

  if(vBuf.size() - iChar < nNbrs*coordSize) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: buffer is too small to extract a Neighborhood" << endl;
    return false;
  }
  pRPL::WeightedCellCoord coord;
  vector<pRPL::WeightedCellCoord> vNbrCoords;
  for(int iNbr = 0; iNbr < nNbrs; iNbr++) {
    memcpy(&coord, &(vBuf[iChar]), coordSize);
    vNbrCoords.push_back(WeightedCellCoord(coord));
    iChar += coordSize;
  }
  return init(vNbrCoords, edgeOption);
}

ostream& pRPL::
operator<<(ostream &os,
           const pRPL::Neighborhood &nbrhd) {
  os << nbrhd.edgeOption() << endl;
  int nNbrs = nbrhd.size();
  os << nNbrs << endl;
  for(int iNbr = 0; iNbr < nNbrs; iNbr++) {
    os << nbrhd[iNbr];
    os << endl;
  }
  return os;
}

istream& pRPL::
operator>>(istream &is,
           pRPL::Neighborhood &nbrhd) {
  pRPL::EdgeOption edgeOption;
  vector<pRPL::WeightedCellCoord> vNbrCoords;
  pRPL::WeightedCellCoord coord;
  int iEdgeOption, nNbrs;
  
  is >> iEdgeOption;
  if(iEdgeOption >= pRPL::FORBID_VIRTUAL_EDGES &&
     iEdgeOption <= pRPL::CUSTOM_VIRTUAL_EDGES) {
    edgeOption = (pRPL::EdgeOption)iEdgeOption;
  }
  else {
    edgeOption = pRPL::FORBID_VIRTUAL_EDGES;
  }

  is >> nNbrs;
  for(int iNbr = 0; iNbr < nNbrs; iNbr++) {
    is >> coord;
    vNbrCoords.push_back(pRPL::WeightedCellCoord(coord));
  }
  nbrhd.init(vNbrCoords, edgeOption);
  
  return is;
}

