#include "prpl-smplDcmp.h"

bool pRPL::SmplDcmp::
_checkDcmpPrmtrs(list<pRPL::SubCellspaceInfo> &lSubspcInfos,
                 const pRPL::SpaceDims &glbDims,
                 int nRowSubspcs, int nColSubspcs,
                 const pRPL::Neighborhood *pNbrhd) const {
  if(!lSubspcInfos.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the vector of SubCellspaceInfo is not empty" \
         << endl;
    return false;
  }
  if(!glbDims.valid()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalide global spatial dimensions (" \
         << glbDims << ")" << endl;
    return false;
  }
  if(nRowSubspcs <= 0 || nRowSubspcs > glbDims.nRows() ||
     nColSubspcs <= 0 || nColSubspcs > glbDims.nCols()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: Invalid numbers of sub-cellspaces (" \
         << nRowSubspcs << ", " << nColSubspcs \
         << "), given the global spatial dimensions (" \
         << glbDims.nRows() << ", " << glbDims.nCols() << ")" \
         << endl;
    return false;
  }
  if(pNbrhd != NULL) {
    if(pNbrhd->isEmpty()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: empty Neighborhood for decomposition" \
           << endl;
      return false;
    }
    else if(pNbrhd->edgeOption() == pRPL::FORBID_VIRTUAL_EDGES &&
            (pNbrhd->nRows() > glbDims.nRows() ||
             pNbrhd->nCols() > glbDims.nCols())) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: the Neighborhood's size (" \
           << pNbrhd->nRows() << ", " << pNbrhd->nCols() \
           << ") is larger than the global spatial dimensions (" \
           << glbDims.nRows() << ", " << glbDims.nCols() << ")" \
           << endl;
      return false;
    }
  }
  return true;
}

void pRPL::SmplDcmp::
_smpl1DDcmp(long &subBegin, long &subEnd,
            long glbBegin, long glbEnd,
            int nSubspcs, int iSubspc) const {
  long glbRange = glbEnd - glbBegin + 1;
  long remainder = glbRange % nSubspcs;
  long blockSize;

  if(remainder == 0) {
    blockSize = glbRange / nSubspcs;
    subBegin = glbBegin + iSubspc * blockSize;
    subEnd = subBegin + blockSize - 1;
  }
  else {
    blockSize = glbRange / nSubspcs + 1;
    if(iSubspc < remainder) {
      subBegin = glbBegin + iSubspc * blockSize;
      subEnd = subBegin + blockSize - 1;
    }
    else {
      if(iSubspc == remainder) {
        subBegin = glbBegin+ iSubspc * blockSize;
        subEnd = subBegin + blockSize - 2;
      }
      else {
        subBegin = glbBegin
                   + remainder * blockSize
                   + (iSubspc - remainder)*(blockSize-1);
        subEnd = subBegin + blockSize - 2;
      }
    }
  }
}

bool pRPL::SmplDcmp::
rowDcmp(list<pRPL::SubCellspaceInfo> &lSubspcInfos,
        const pRPL::SpaceDims &glbDims,
        int nSubspcs,
        const pRPL::Neighborhood *pNbrhd) const {
  if(!_checkDcmpPrmtrs(lSubspcInfos, glbDims,
                       nSubspcs, 1, pNbrhd)) {
    return false;
  }

  lSubspcInfos.clear();

  pRPL::CoordBR glbWorkBR;
  if(pNbrhd != NULL) {
    if(!pNbrhd->calcWorkBR(glbWorkBR, glbDims)) {
      return false;
    }
  }
  else {
    glbWorkBR.nwCorner(0, 0);
    glbWorkBR.seCorner(glbDims.nRows()-1, glbDims.nCols()-1);
  }
  if(nSubspcs > glbWorkBR.nRows()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: too many sub-cellspaces (" \
         << nSubspcs << "), given the global working dimensions (" \
         << glbWorkBR.nRows() << ", " << glbWorkBR.nCols() << ")" \
         << endl;
    return false;
  }

  pRPL::DomDcmpType dcmpType = pRPL::ROW_DCMP;
  long glbBegin = glbWorkBR.nwCorner().iRow();
  long glbEnd = glbWorkBR.seCorner().iRow();

  for(int iSubspc = 0; iSubspc < nSubspcs; iSubspc++) {
    long subBegin, subEnd;
    _smpl1DDcmp(subBegin, subEnd,
                glbBegin, glbEnd,
                nSubspcs, iSubspc);

    pRPL::CellCoord nwCorner, seCorner;
    if(pNbrhd != NULL) {
      nwCorner.iRow(subBegin + pNbrhd->minIRow());
      seCorner.iRow(subEnd + pNbrhd->maxIRow());
      if(pNbrhd->edgeOption() == pRPL::NODATA_VIRTUAL_EDGES ||
         pNbrhd->edgeOption() == pRPL::CUSTOM_VIRTUAL_EDGES) {
        if(iSubspc == 0) {
          nwCorner.iRow(subBegin);
        }
        if(iSubspc == nSubspcs - 1) {
          seCorner.iRow(subEnd);
        }
      }
    }
    else {
      nwCorner.iRow(subBegin);
      seCorner.iRow(subEnd);
    }
    nwCorner.iCol(0);
    seCorner.iCol(glbDims.nCols() - 1);

    pRPL::CoordBR subMBR(nwCorner, seCorner);
    pRPL::SpaceDims subDims(subMBR.nRows(), subMBR.nCols());

    pRPL::CoordBR subWorkBR;
    if(pNbrhd != NULL) {
      if(!pNbrhd->calcWorkBR(subWorkBR, subDims)) {
        return false;
      }
      if(pNbrhd->edgeOption() == pRPL::NODATA_VIRTUAL_EDGES ||
         pNbrhd->edgeOption() == pRPL::CUSTOM_VIRTUAL_EDGES) {
        if(iSubspc != 0 && iSubspc != nSubspcs - 1) {
          subWorkBR.nwCorner(-(pNbrhd->minIRow()),
                             subWorkBR.nwCorner().iCol());
          subWorkBR.seCorner(subDims.nRows() - pNbrhd->maxIRow() - 1,
                             subWorkBR.seCorner().iCol());
        }
        if(iSubspc == 0) {
          subWorkBR.seCorner(subDims.nRows() - pNbrhd->maxIRow() - 1,
                             subWorkBR.seCorner().iCol());
        }
        if(iSubspc == nSubspcs - 1) {
          subWorkBR.nwCorner(-(pNbrhd->minIRow()),
                             subWorkBR.nwCorner().iCol());
        }
      }
    }
    else {
      subWorkBR.nwCorner(0, 0);
      subWorkBR.seCorner(subDims.nRows()-1, subDims.nCols()-1);
    }

    vector<pRPL::IntVect> mNbrSpcIDs(2);
    if(iSubspc != 0) {
      mNbrSpcIDs[pRPL::UPPER_DIR].push_back(iSubspc - 1);
    }
    if(iSubspc != nSubspcs-1) {
      mNbrSpcIDs[pRPL::LOWER_DIR].push_back(iSubspc + 1);
    }

    lSubspcInfos.push_back(pRPL::SubCellspaceInfo(iSubspc, dcmpType, glbDims,
                                                  subMBR, subWorkBR, mNbrSpcIDs));

    /*
    cout << iSubspc << " glbDims[" << glbDims << "]; subMBR[" << subMBR \
               << "]; subWorkBR[" << subWorkBR << "]" << endl;
    */

  } // End of for(iSubspc)

  return true;
}

bool pRPL::SmplDcmp::
colDcmp(list<pRPL::SubCellspaceInfo> &lSubspcInfos,
        const pRPL::SpaceDims &glbDims,
        int nSubspcs,
        const pRPL::Neighborhood *pNbrhd) const {
  if(!_checkDcmpPrmtrs(lSubspcInfos, glbDims,
                       1, nSubspcs, pNbrhd)) {
    return false;
  }

  lSubspcInfos.clear();

  pRPL::CoordBR glbWorkBR;
  if(pNbrhd != NULL) {
    if(!pNbrhd->calcWorkBR(glbWorkBR, glbDims)) {
      return false;
    }
  }
  else {
    glbWorkBR.nwCorner(0, 0);
    glbWorkBR.seCorner(glbDims.nRows()-1, glbDims.nCols()-1);
  }
  if(nSubspcs > glbWorkBR.nCols()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: too many sub-cellspaces (" \
         << nSubspcs << "), given the global working dimensions (" \
         << glbWorkBR.nRows() << ", " << glbWorkBR.nCols() << ")" \
         << endl;
    return false;
  }

  pRPL::DomDcmpType dcmpType = pRPL::COL_DCMP;
  long glbBegin = glbWorkBR.nwCorner().iCol();
  long glbEnd = glbWorkBR.seCorner().iCol();

  for(int iSubspc = 0; iSubspc < nSubspcs; iSubspc++) {
    long subBegin, subEnd;
    _smpl1DDcmp(subBegin, subEnd,
                glbBegin, glbEnd,
                nSubspcs, iSubspc);

    pRPL::CellCoord nwCorner, seCorner;
    if(pNbrhd != NULL) {
      nwCorner.iCol(subBegin + pNbrhd->minICol());
      seCorner.iCol(subEnd + pNbrhd->maxICol());
      if(pNbrhd->edgeOption() == pRPL::NODATA_VIRTUAL_EDGES ||
         pNbrhd->edgeOption() == pRPL::CUSTOM_VIRTUAL_EDGES) {
        if(iSubspc == 0) {
          nwCorner.iCol(subBegin);
        }
        if(iSubspc == nSubspcs - 1) {
          seCorner.iCol(subEnd);
        }
      }
    }
    else {
      nwCorner.iCol(subBegin);
      seCorner.iCol(subEnd);
    }
    nwCorner.iRow(0);
    seCorner.iRow(glbDims.nRows() - 1);

    pRPL::CoordBR subMBR(nwCorner, seCorner);
    pRPL::SpaceDims subDims(subMBR.nRows(), subMBR.nCols());

    pRPL::CoordBR subWorkBR;
    if(pNbrhd != NULL) {
      if(!pNbrhd->calcWorkBR(subWorkBR, subDims)) {
        return false;
      }
      if(pNbrhd->edgeOption() == pRPL::NODATA_VIRTUAL_EDGES ||
         pNbrhd->edgeOption() == pRPL::CUSTOM_VIRTUAL_EDGES) {
        if(iSubspc != 0 && iSubspc != nSubspcs - 1) {
          subWorkBR.nwCorner(subWorkBR.nwCorner().iRow(),
                             -(pNbrhd->minICol()));
          subWorkBR.seCorner(subWorkBR.seCorner().iRow(),
                             subDims.nCols() - pNbrhd->maxICol() - 1);
        }
        if(iSubspc == 0) {
          subWorkBR.seCorner(subWorkBR.seCorner().iRow(),
                             subDims.nCols() - pNbrhd->maxICol() - 1);
        }
        if(iSubspc == nSubspcs - 1) {
          subWorkBR.nwCorner(subWorkBR.nwCorner().iRow(),
                             -(pNbrhd->minICol()));
        }
      }
    }
    else {
      subWorkBR.nwCorner(0, 0);
      subWorkBR.seCorner(subDims.nRows()-1, subDims.nCols()-1);
    }

    vector<pRPL::IntVect> mNbrSpcIDs(2);
    if(iSubspc != 0) {
      mNbrSpcIDs[pRPL::LEFT_DIR].push_back(iSubspc - 1);
    }
    if(iSubspc != nSubspcs-1) {
      mNbrSpcIDs[pRPL::RIGHT_DIR].push_back(iSubspc + 1);
    }

    lSubspcInfos.push_back(pRPL::SubCellspaceInfo(iSubspc, dcmpType, glbDims,
                                                  subMBR, subWorkBR, mNbrSpcIDs));
  } // End of for(iSubspc)

  return true;
}

bool pRPL::SmplDcmp::
blkDcmp(list<pRPL::SubCellspaceInfo> &lSubspcInfos,
        const pRPL::SpaceDims &glbDims,
        int nRowSubspcs,
        int nColSubspcs,
        const pRPL::Neighborhood *pNbrhd) const {
  if(!_checkDcmpPrmtrs(lSubspcInfos, glbDims, nRowSubspcs, nColSubspcs, pNbrhd)) {
      return false;
    }

  lSubspcInfos.clear();

  pRPL::CoordBR glbWorkBR;
  if(pNbrhd != NULL) {
    if(!pNbrhd->calcWorkBR(glbWorkBR, glbDims)) {
      return false;
    }
  } // end of if(pNbrhd != NULL)
  else {
    glbWorkBR.nwCorner(0, 0);
    glbWorkBR.seCorner(glbDims.nRows()-1, glbDims.nCols()-1);
  }
  if(nRowSubspcs > glbWorkBR.nRows() ||
     nColSubspcs > glbWorkBR.nCols()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: too many sub-cellspaces (" \
         << nRowSubspcs << ", " << nColSubspcs \
         << "), given the global working dimensions (" \
         << glbWorkBR.nRows() << ", " << glbWorkBR.nCols() << ")" \
         << endl;
    return false;
  }

  pRPL::DomDcmpType dcmpType = pRPL::BLK_DCMP;
  long glbRowBegin = glbWorkBR.nwCorner().iRow();
  long glbRowEnd = glbWorkBR.seCorner().iRow();
  long glbColBegin = glbWorkBR.nwCorner().iCol();
  long glbColEnd = glbWorkBR.seCorner().iCol();

  for(int iRowSubspc = 0; iRowSubspc < nRowSubspcs; iRowSubspc++) {
    long subRowBegin, subRowEnd;
    _smpl1DDcmp(subRowBegin, subRowEnd,
                glbRowBegin, glbRowEnd,
                nRowSubspcs, iRowSubspc);

    for(int iColSubspc = 0; iColSubspc < nColSubspcs; iColSubspc++) {
      long subColBegin, subColEnd;
      _smpl1DDcmp(subColBegin, subColEnd,
                  glbColBegin, glbColEnd,
                  nColSubspcs, iColSubspc);

      pRPL::CellCoord nwCorner, seCorner;
      if(pNbrhd != NULL) {
        nwCorner.iRow(subRowBegin + pNbrhd->minIRow());
        nwCorner.iCol(subColBegin + pNbrhd->minICol());
        seCorner.iRow(subRowEnd + pNbrhd->maxIRow());
        seCorner.iCol(subColEnd + pNbrhd->maxICol());
        if(pNbrhd->edgeOption() == pRPL::NODATA_VIRTUAL_EDGES ||
           pNbrhd->edgeOption() == pRPL::CUSTOM_VIRTUAL_EDGES) {
          if(iRowSubspc == 0) {
            nwCorner.iRow(subRowBegin);
          }
          if(iRowSubspc == nRowSubspcs - 1) {
            seCorner.iRow(subRowEnd);
          }
          if(iColSubspc == 0) {
            nwCorner.iCol(subColBegin);
          }
          if(iColSubspc == nColSubspcs - 1) {
            seCorner.iCol(subColEnd);
          }
        }
      }
      else {
        nwCorner.iRow(subRowBegin);
        nwCorner.iCol(subColBegin);
        seCorner.iRow(subRowEnd);
        seCorner.iCol(subColEnd);
      }

      pRPL::CoordBR subMBR(nwCorner, seCorner);
      pRPL::SpaceDims subDims(subMBR.nRows(), subMBR.nCols());

      pRPL::CoordBR subWorkBR;
      if(pNbrhd != NULL) {
        if(!pNbrhd->calcWorkBR(subWorkBR, subDims)) {
          return false;
        }
        if(pNbrhd->edgeOption() == pRPL::NODATA_VIRTUAL_EDGES ||
           pNbrhd->edgeOption() == pRPL::CUSTOM_VIRTUAL_EDGES) {
          if(iRowSubspc != 0 && iRowSubspc != nRowSubspcs - 1) {
            subWorkBR.nwCorner(-(pNbrhd->minIRow()),
                               subWorkBR.nwCorner().iCol());
            subWorkBR.seCorner(subDims.nRows() - pNbrhd->maxIRow() - 1,
                               subWorkBR.seCorner().iCol());
          }
          if(iRowSubspc == 0) {
            subWorkBR.seCorner(subDims.nRows() - pNbrhd->maxIRow() - 1,
                               subWorkBR.seCorner().iCol());
          }
          if(iRowSubspc == nRowSubspcs - 1) {
            subWorkBR.nwCorner(-(pNbrhd->minIRow()),
                               subWorkBR.nwCorner().iCol());
          }
          if(iColSubspc != 0 && iColSubspc != nColSubspcs - 1) {
            subWorkBR.nwCorner(subWorkBR.nwCorner().iRow(),
                               -(pNbrhd->minICol()));
            subWorkBR.seCorner(subWorkBR.seCorner().iRow(),
                               subDims.nCols() - pNbrhd->maxICol() - 1);
          }
          if(iColSubspc == 0) {
            subWorkBR.seCorner(subWorkBR.seCorner().iRow(),
                               subDims.nCols() - pNbrhd->maxICol() - 1);
          }
          if(iColSubspc == nColSubspcs - 1) {
            subWorkBR.nwCorner(subWorkBR.nwCorner().iRow(),
                               -(pNbrhd->minICol()));
          }
        }
      }
      else {
        subWorkBR.nwCorner(0, 0);
        subWorkBR.seCorner(subDims.nRows()-1, subDims.nCols()-1);
      }

      vector<pRPL::IntVect> mNbrSpcIDs(8);
      for(int iDir = 0; iDir < 8; iDir++) {
        int iRowNbrSpc = iRowSubspc + MooreNbrLocs[2*iDir];
        int iColNbrSpc = iColSubspc + MooreNbrLocs[2*iDir + 1];
        if(iRowNbrSpc >= 0 && iRowNbrSpc < nRowSubspcs &&
           iColNbrSpc >= 0 && iColNbrSpc < nColSubspcs) {
          int iNbrSpc = iRowNbrSpc * nColSubspcs + iColNbrSpc;
          mNbrSpcIDs[iDir].push_back(iNbrSpc);
        }
      }

      int iSubspc = iRowSubspc * nColSubspcs + iColSubspc;
      lSubspcInfos.push_back(pRPL::SubCellspaceInfo(iSubspc, dcmpType, glbDims,
                                                    subMBR, subWorkBR, mNbrSpcIDs));

    } // end of for(iColSubspc)
  } // end of for(iRowSubspc)

  return true;
}
