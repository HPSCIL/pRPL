#ifndef PRPL_NEIGHBORHOOD_H
#define PRPL_NEIGHBORHOOD_H

/***************************************************************************
* prpl-neighborhood.h
*
* Project: pRPL, v 2.0
* Purpose: Header file for class pRPL::Neighborhood
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

namespace pRPL {
  static const int MooreNbrLocs[16] = {-1, 0,
                                       -1, 1,
                                       0, 1,
                                       1, 1,
                                       1, 0,
                                       1, -1,
                                       0, -1,
                                       -1, -1};

  class Neighborhood {
    public:
      /* Constructor and destructor */
      Neighborhood(const char *aName = NULL,
                   pRPL::EdgeOption edgeOption = FORBID_VIRTUAL_EDGES,
                   double virtualEdgeVal = 0.0);
      Neighborhood(const pRPL::Neighborhood &rhs);
      ~Neighborhood();
      
      /* Initiate */
      bool init(const vector<pRPL::CellCoord> &vNbrCoords,
                double weight = 1.0,
                pRPL::EdgeOption edgeOption = FORBID_VIRTUAL_EDGES,
                double virtualEdgeVal = 0.0);
      bool init(const vector<pRPL::WeightedCellCoord> &vNbrCoords,
                pRPL::EdgeOption edgeOption = FORBID_VIRTUAL_EDGES,
                double virtualEdgeVal = 0.0);
      bool initSingleCell(double weight = 1.0);
      bool initVonNeumann(double weight = 1.0,
                          pRPL::EdgeOption edgeOption = FORBID_VIRTUAL_EDGES,
                          double virtualEdgeVal = 0.0);
      bool initMoore(long nEdgeLength = 3,
                     double weight = 1.0,
                     pRPL::EdgeOption edgeOption = FORBID_VIRTUAL_EDGES,
                     double virtualEdgeVal = 0.0);
      void clear();
      
      /* Operators */
      pRPL::WeightedCellCoord& operator[](int iNbr);
      const pRPL::WeightedCellCoord& operator[](int iNbr) const;
      pRPL::Neighborhood& operator=(const pRPL::Neighborhood &rhs);

      /* Information */
      const char* name() const;
      void name(const char *aName);
      bool isEmpty() const;
      int size() const;
      bool isEquallyWeighted(double &weight) const;
      long minIRow() const;
      long minICol() const;
      long maxIRow() const;
      long maxICol() const;
      long nRows() const;
      long nCols() const;
      const pRPL::CoordBR* getMBR() const;
      pRPL::EdgeOption edgeOption() const;
      const double* virtualEdgeVal() const;
      bool virtualEdgeVal(const double veVal);
      
      /* Neighbor accessing */
      bool hasNbrs(pRPL::MeshDir dir) const;
      const pRPL::IntVect* nbrIDs(pRPL::MeshDir dir) const;
      
      /* Cellspace-related methods */
      bool calcWorkBR(pRPL::CoordBR &workBR,
                      const pRPL::SpaceDims &dims) const;

      /* Compression */
      void add2Buf(pRPL::CharVect &vBuf) const;
      bool fromBuf(const pRPL::CharVect &vBuf,
                   int &iChar);
      
    private:
      bool _checkNbrCoords(const vector<pRPL::CellCoord> &vNbrCoords) const;
      bool _validNumDirects() const;

    private:
      string _name;
      vector<pRPL::WeightedCellCoord> _vNbrs;
      pRPL::CoordBR _MBR;
      vector<pRPL::IntVect> _mNbrIDMap;
      pRPL::EdgeOption _edgeOption;
      double *_pVirtualEdgeVal;
  };

  ostream& operator<<(ostream &os,
                      const pRPL::Neighborhood &nbrhd);
  
  istream& operator>>(istream &is,
                      pRPL::Neighborhood &nbrhd);
};

#endif
