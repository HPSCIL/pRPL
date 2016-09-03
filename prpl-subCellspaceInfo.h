#ifndef PRPL_SUBCELLSPACEINFO_H
#define PRPL_SUBCELLSPACEINFO_H

/***************************************************************************
* prpl-subCellspaceInfo.h
*
* Project: pRPL, v 2.2
* Purpose: Header file for class pRPL::SubCellspaceInfo
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

#include "prpl-basicTypes.h"
#include "prpl-neighborhood.h"

namespace pRPL {
  class SubCellspaceInfo {
    public:
      SubCellspaceInfo();
      SubCellspaceInfo(int id,
                       pRPL::DomDcmpType domDcmpType,
                       const pRPL::SpaceDims &glbDims,
                       const pRPL::CoordBR &MBR,
                       const pRPL::CoordBR &workBR,
                       const vector<pRPL::IntVect> &mNbrSpcIDs);
      SubCellspaceInfo(const pRPL::SubCellspaceInfo &rhs);
      SubCellspaceInfo(const pRPL::IntVect &vInfoPack,
                       pRPL::IntVctItr &iVal);

      ~SubCellspaceInfo() {}

      pRPL::SubCellspaceInfo& operator=(const pRPL::SubCellspaceInfo &rhs);
      bool operator==(const pRPL::SubCellspaceInfo &rhs) const;
      bool operator!=(const pRPL::SubCellspaceInfo &rhs) const;

      int id() const;
      pRPL::DomDcmpType domDcmpType() const;
      long nGlbRows() const;
      long nGlbCols() const;
      const pRPL::SpaceDims& glbDims() const;
      long nRows() const;
      long nCols() const;
      const pRPL::SpaceDims& dims() const;
      long iRowBegin() const;
      long iColBegin() const;
      long iRowEnd() const;
      long iColEnd() const;
      const pRPL::CoordBR& MBR() const;
      long iRowWorkBegin() const;
      long iColWorkBegin() const;
      long iRowWorkEnd() const;
      long iColWorkEnd() const;
      const pRPL::CoordBR& workBR() const;
      double sizeRatio() const;
      
      bool validGlbIdx(long glbIdx,
                       bool warning = true) const;
      bool validGlbCoord(const pRPL::CellCoord &glbCoord,
                         bool warning = true) const;
      bool validGlbCoord(long iRowGlb, long iColGlb,
                         bool warning = true) const;
      const pRPL::CellCoord glbIdx2glbCoord(long glbIdx) const;
      long glbCoord2glbIdx(const pRPL::CellCoord &glbCoord) const;
      long glbCoord2glbIdx(long iRowGlb, long iColGlb) const;

      bool validIdx(long idx,
                    bool warning = true) const;
      bool validCoord(const pRPL::CellCoord &coord,
                      bool warning = true) const;
      bool validCoord(long iRow, long iCol,
                      bool warning = true) const;
      long coord2idx(const pRPL::CellCoord &coord) const;
      long coord2idx(long iRow, long iCol) const;
      const pRPL::CellCoord idx2coord(long idx) const;

      const pRPL::CellCoord lclCoord2glbCoord(const pRPL::CellCoord& lclCoord) const;
      const pRPL::CellCoord lclCoord2glbCoord(long iRowLcl, long iColLcl) const;
      const pRPL::CellCoord glbCoord2lclCoord(const pRPL::CellCoord& glbCoord) const;
      const pRPL::CellCoord glbCoord2lclCoord(long iRowGlb, long iColGlb) const;
      long lclCoord2glbIdx(const pRPL::CellCoord &lclCoord) const;
      long lclCoord2glbIdx(long iRowLcl, long iColLcl) const;
      const pRPL::CellCoord glbIdx2lclCoord(long glbIdx) const;
      long glbIdx2lclIdx(long glbIdx) const;
      long lclIdx2glbIdx(long lclIdx) const;

      int nNbrDirs() const;
      int nEdges() const;
      int nTotNbrs() const;
      bool validSpcDir(int iDir) const;
      bool hasNbrs(int iDir) const;
      int nNbrs(int iDir) const;
      const pRPL::IntVect& nbrSubSpcIDs(int iDir) const;
      int nbrDir(int nbrID) const;
      pRPL::MeshDir spcDir2MeshDir(int iDir) const;
      int oppositeDir(int iDir) const;

      bool calcBRs(bool onlyUpdtCtrCell = true,
                   const pRPL::Neighborhood *pNbrhd = NULL,
                   const list<pRPL::SubCellspaceInfo> *plSubspcInfos = NULL);
      const pRPL::CoordBR* interiorBR() const;
      const pRPL::CoordBR* edgeBR(int iDir) const;
      const pRPL::CoordBR* sendBR(int iDir,
                                  int iNbr = 0) const;

      /*
      void add2IntVect(pRPL::IntVect &vInfoPack) const;
      bool fromIntVect(const pRPL::IntVect &vInfoPack,
                       pRPL::IntVctItr &iVal);
      */

    protected:
      void _initBRs();
      bool _rowwiseBRs(bool onlyUpdtCtrCell = true,
                       const pRPL::Neighborhood *pNbrhd = NULL);
      bool _colwiseBRs(bool onlyUpdtCtrCell = true,
                       const pRPL::Neighborhood *pNbrhd = NULL);
      bool _blockwiseBRs(bool onlyUpdtCtrCell = true,
                         const pRPL::Neighborhood *pNbrhd = NULL,
                         const list<pRPL::SubCellspaceInfo> *plSubspcInfos = NULL);
      const pRPL::SubCellspaceInfo* _findSubspcInfo(const list<pRPL::SubCellspaceInfo> *plSubspcInfos,
                                                    int subspcGlbID) const;

    protected:
      int _id;
      pRPL::DomDcmpType _domDcmpType;
      pRPL::SpaceDims _glbDims;
      pRPL::CoordBR _MBR; /* MBR is in global coordinates */
      pRPL::SpaceDims _dims;
      pRPL::CoordBR _workBR; /* workBR is in local coordinates */
      vector<pRPL::IntVect> _mNbrSpcIDs;

      vector<vector<CoordBR> > _mSendBRs;
      vector<CoordBR> _vEdgeBRs;
      CoordBR _interiorBR;
  };
};

#endif
