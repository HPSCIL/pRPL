#ifndef PRPL_BASICTYPES_H
#define PRPL_BASICTYPES_H

/***************************************************************************
* prpl-basicTypes.h
*
* Project: pRPL, v 2.0
* Purpose: Header file for class pRPL::DomDcmpType, pRPL::DomDcmpMethod,
*          pRPL::SglColDir, pRPL::SglRowDir, pRPL::MeshDir, pRPL::PrimeDir,
*          pRPL::SpaceDims, pRPL::CellCoord, and pRPL::CoordBR
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <typeinfo>
#include <utility>
#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cwchar>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "pgtiol-basicTypes.h"
#include "pgtiol-dataset.h"

using namespace std;

namespace pRPL {
  /***********************************************
   * Basic Data Types                            *
   ***********************************************/  
  void gdalType2cppType(string &cppTypeName,
                        size_t &cppTypeSize,
                        GDALDataType gdalType);

  void cppType2gdalType(GDALDataType &gdalType,
                        const string &cppTypeName);

  bool initDefVal(void *&pVal,
                  const string &cppTypeName);

  /***********************************************
   * Neighborhood Virtual Edge Options           *
   * FORBID_VIRTUAL_EDGES: no virtual edges are  *
   *                       allowed               *
   * NODATA_VIRTUAL_EDGES: assign NoData value   *
   *                       to virtual edge cells *
   ***********************************************/
  enum EdgeOption{FORBID_VIRTUAL_EDGES = 0,
                  NODATA_VIRTUAL_EDGES,
                  CUSTOM_VIRTUAL_EDGES};

  /***********************************************
   * Transition "Evaluate" Returns               *
   ***********************************************/  
  enum EvaluateReturn{EVAL_FAILED = 0,
  	                  EVAL_SUCCEEDED,
  	                  EVAL_TERMINATED};
  
  /***********************************************
   * Transition "Evaluate" Types                 *
   ***********************************************/
  enum EvaluateType{EVAL_NONE = 0,
                    EVAL_ALL,
                    EVAL_RANDOMLY,
                    EVAL_SELECTED};

  /***********************************************
   * Transition "Evaluate" BRs for a SubCellspace*
   ***********************************************/
  enum EvaluateBR{EVAL_WORKBR = 0,
                  EVAL_EDGES,
                  EVAL_INTERIOR};

  /***********************************************
  * Special Process Group IDs                    *
  ***********************************************/
  enum PrcGrpID{ROOT_GRP = -1,
                MASTER_GRP = -99};
  
  /***********************************************
  * Domain Decomposition Types                   *
  * NOTE: Quad-tree Decomposition is a special   *
  *       case of Block-wise Decomposition       *
  ***********************************************/
  enum DomDcmpType{NON_DCMP = 0,
                   ROW_DCMP,
                   COL_DCMP,
                   BLK_DCMP};

  /***********************************************
  * Domain Decomposition Methods                 *
  ***********************************************/
  enum DomDcmpMethod{SNGL_PRC = 0,
                     SMPL_ROW,
                     SMPL_COL,
                     SMPL_BLK,
                     WKLD_ROW,
                     WKLD_COL,
                     WKLD_BLK,
                     WKLD_QTB};

  /***********************************************
  * Sub-cellspace-to-Process Mapping Methods     *
  * CYLC_MAP: Cyclic Mapping;                     *
  * BKFR_MAP: Back-n-Forth Mapping               *
  ***********************************************/
  enum MappingMethod{CYLC_MAP = 0,
                     BKFR_MAP};

  /***********************************************
  * Output-Stream Operators for                  *
  * Decomposition Method and Mapping Method      *
  ***********************************************/
  ostream& operator<<(ostream &os, const DomDcmpMethod &dcmpMethod);
  ostream& operator<<(ostream &os, const MappingMethod &mapMethod);
  
  /***********************************************
  * Load-balancing Methods                       *
  * STC_LDBLC: Static Load-balancing;            *
  * DYN_LDBLC: Dynamic Load-balancing            *
  ***********************************************/
  enum LdBlcMethod{STC_LDBLC = 0,
                   DYN_LDBLC};

  /***********************************************
  * Task Signals                              *
  ***********************************************/
  enum TaskSignal{NEWTASK_SGNL = -999,
                  SUBMIT_SGNL = -9999,
                  QUIT_SGNL = -99999};

  /***********************************************
  * Message Tags                                 *
  * DSTRBT_TAG:   Distribute                     *
  * TRSFER_TAG:   Transfer                       *
  * REQUEST_TAG:  Request                        *
  * SIZE_TAG:     Size                           *
  * EXCHANGE_TAG: Exchange                       *
  ***********************************************/
  enum MsgTag{DSTRBT_TAG = 999,
              SUBMIT_TAG = 9999,
              REQUEST_TAG = 99999,
              SIZE_TAG = 999999,
              EXCHANGE_TAG = 9999999};

  /***********************************************
   * Transfer Data Type                          *
   * SPACE_TRANS:   Transfer a (Sub)Cellspace    *
   * CELL_TRANS:    Transfer Cells               *
   ***********************************************/
  enum TransferDataType{SPACE_TRANSFDATA = 0,
                        CELL_TRANSFDATA};

  /***********************************************
   * Transfer Type                               *
   * BCST_TRANSF:    Broadcasting Transfer       *
   * SEND_TRANSF:    Sending Transfer            *
   * RECV_TRANSF:    Receiving Transfer          *
   ***********************************************/
  enum TransferType{BCST_TRANSF = 0,
                    SEND_TRANSF,
                    RECV_TRANSF};

  /***************************************************
   * Reading Options                                 *
   * PARA_READING:    Parallel Reading               *
   * CENT_READING:    Centralized Reading            *
   * CENTDEL_READING: Centralized Reading &          *
   *                  Delete SubCellspaces on Master *
   * PGT_READING:     Parallel GeoTiff I/O Library   *
   ***************************************************/
  enum ReadingOption{NO_READING = -1,
                     PARA_READING = 0,
                     CENT_READING,
                     CENTDEL_READING,
                     PGT_READING};

  /***************************************************
   * Writing Options                                 *
   * NO_WRITING:      No Writing                     *
   * PARA_WRITING:    Parallel Writing               *
   * PARADEL_WRITING: Parallel Writing &             *
   *                  Delete SubCellspaces on Workers*
   * CENT_WRITING:    Centralized Writing            *
   * CENTDEL_WRITING: Centralized Writing &          *
   *                  Delete SubCellspaces on Master *
   *                  and Workers                    *
   * PGT_WRITING:     Parallel GeoTiff I/O Library   *
   * PGTDEL_WRITING:  Parallel GeoTiff I/O Library & *
   *                  Delete SubCellspaces on Workers*
   ***************************************************/
  enum WritingOption{NO_WRITING = -1,
                     PARA_WRITING = 0,
                     PARADEL_WRITING,
                     CENT_WRITING,
                     CENTDEL_WRITING,
                     PGT_WRITING,
                     PGTDEL_WRITING};

  /***********************************************
  * Vertical Directions                          *
  ***********************************************/
  enum VerticalDir{UPPER_DIR = 0,
                   LOWER_DIR};

  /***********************************************
  * Horizontal Directions                        *
  ***********************************************/
  enum HorizontalDir{LEFT_DIR = 0,
                     RIGHT_DIR};

  /***********************************************
  * Directions in a Mesh (8 directions)          *
  ***********************************************/
  enum MeshDir{NONE_DIR = -1,
               NORTH_DIR = 0,
               NORTHEAST_DIR,
               EAST_DIR,
               SOUTHEAST_DIR,
               SOUTH_DIR,
               SOUTHWEST_DIR,
               WEST_DIR,
               NORTHWEST_DIR,
               SPECIAL1_DIR,
               SPECIAL2_DIR};
  
  /***********************************************
  * Function to return the opposite direction    *
  * of a given mesh direction                    *
  ***********************************************/
  pRPL::MeshDir oppositeDir(pRPL::MeshDir dir);
  
  /***********************************************
  * Prime Directions in a Mesh (4 directions)    *
  ***********************************************/
  enum PrimeDir{NONE_PDIR = -1,
                NORTH_PDIR = 0,
                EAST_PDIR,
                SOUTH_PDIR,
                WEST_PDIR};

  /***********************************************
   * Function to return the opposite direction    *
   * of a given prime direction                   *
   ***********************************************/
  pRPL::PrimeDir oppositePrmDir(pRPL::PrimeDir dir);

  /***********************************************
  * Neighborhood types                           *
  ***********************************************/
  enum BasicNeighborhoodType{SINGLE_CELL_NBRHD = 0,
                             VON_NEUMANN_NBRHD,
                             MOORE_NBRHD};

  ostream& operator<<(ostream &os, const pRPL::BasicNeighborhoodType &val);
  istream& operator>>(istream &is, pRPL::BasicNeighborhoodType &val);

  /***********************************************
  * Vector of Integer Members                    *
  ***********************************************/
  typedef vector<int> IntVect;
  typedef vector<int>::iterator IntVctItr;
  ostream& operator<<(ostream &os, const pRPL::IntVect &vInt);

  /***********************************************
   * Vector of Long Integer Members                    *
   ***********************************************/
  typedef vector<long> LongVect;
  typedef vector<long>::iterator LongVctItr;
  ostream& operator<<(ostream &os, const pRPL::LongVect &vLong);
  
  /***********************************************
  * Vector of Double-floating Members            *
  ***********************************************/
  typedef vector<double> DblVect;
  typedef vector<double>::iterator DblVectItr;
  ostream& operator<<(ostream &os, const pRPL::DblVect &vDbl);

  /***********************************************
  * Vector of Char Members                       *
  ***********************************************/
  typedef vector<char> CharVect;
  typedef vector<char>::iterator CharVctItr;
  ostream& operator<<(ostream &os, const pRPL::CharVect &vChar);

  /***********************************************
  * Integer-Integer Pair                         *
  ***********************************************/
  typedef pair<int, int> IntPair;
  ostream& operator<<(ostream &os, const pRPL::IntPair &pair);

  /***********************************************
  * Integer-Integer Map                          *
  * Integer Key and Integer Value                *
  ***********************************************/
  typedef map<int, int> IntMap;
  typedef map<int, int>::value_type IntMapValType;
  typedef map<int, int>::iterator IntMapItr;
  ostream& operator<<(ostream &os, const pRPL::IntMap &mInt);

  /***********************************************
  * Constants                                    *
  ***********************************************/
  static const int              ERROR_ID = -1;
  static const long             ERROR_COORD = -9999;
  static const int              ERROR_VAL = -9999;
  static const double           EPSINON = 1e-06;

  static const bool             DEFAULT_NODATA_BOOL = false;
  static const unsigned char    DEFAULT_NODATA_UCHAR = 255;
  static const char             DEFAULT_NODATA_CHAR = -128;
  static const unsigned short   DEFAULT_NODATA_USHORT = 65535;
  static const short            DEFAULT_NODATA_SHORT = -32768;
  static const unsigned int     DEFAULT_NODATA_UINT = 4294967295;
  static const int              DEFAULT_NODATA_INT = -2147483648;
  static const unsigned long    DEFAULT_NODATA_ULONG = 4294967295;
  static const long             DEFAULT_NODATA_LONG = -2147483648;
  static const float            DEFAULT_NODATA_FLOAT = -2147483648;
  static const double           DEFAULT_NODATA_DOUBLE = -2147483648;
  static const long double      DEFAULT_NODATA_LDOUBLE = -2147483648;

  /***********************************************
  * SpaceDims class                              *
  * Defining the spatial dimensions of a         *
  * cellspace: Number of Rows; Number of Columns *
  ***********************************************/
  class SpaceDims {
    public:
      SpaceDims()
        :_nRows(0),
         _nCols(0),
         _size(0) {}
      SpaceDims(long nRows, long nCols)
        :_nRows(nRows),
         _nCols(nCols),
         _size(nRows*nCols) {}
      SpaceDims(const pRPL::SpaceDims &rhs)
        :_nRows(rhs._nRows),
         _nCols(rhs._nCols),
         _size(rhs._size) {}
      ~SpaceDims() {}

      pRPL::SpaceDims& operator=(const pRPL::SpaceDims &rhs) {
        if(this != &rhs) {
          _nRows = rhs._nRows;
          _nCols = rhs._nCols;
          _size = rhs._size;
        }
        return *this;
      }
      bool operator==(const pRPL::SpaceDims &rhs) const {
        return (_nRows == rhs._nRows &&
                _nCols == rhs._nCols);
      }
      bool operator!=(const pRPL::SpaceDims &rhs) const {
        return (_nRows != rhs._nRows ||
                _nCols != rhs._nCols);
      }
      void nRows(long numRows) {
        _nRows = numRows;
        _size = _nRows * _nCols;
      }
      void nCols(long numCols) {
        _nCols = numCols;
        _size = _nRows * _nCols;
      }
      long nRows() const {
        return _nRows;
      }
      long nCols() const {
        return _nCols;
      }
      long size() const {
        return _size;
      }
      bool valid() const {
        return (_nRows >= 0 && _nCols >= 0);
      }
      bool isNone() const {
        return (_size <= 0);
      }
      bool validIdx(long idx) const {
        return (idx >= 0 && idx < _size);
      }

      void add2Buf(vector<char> &buf) const {
        size_t bufPtr = buf.size();
        buf.resize(bufPtr + 2*sizeof(long));

        memcpy((void *)&(buf[bufPtr]), &_nCols, sizeof(long));
        bufPtr += sizeof(long);
        memcpy((void *)&(buf[bufPtr]), &_nRows, sizeof(long));
      }

      bool initFromBuf(vector<char> &buf) {
        _nRows = 0;
        _nCols = 0;
        _size = 0;

        int bufPtr = buf.size() - sizeof(long);
        if(bufPtr < 0) {
          cerr << __FILE__ << " function:" << __FUNCTION__ \
              << " Error: the buffer has insufficient memory to extract the number of rows" \
              << endl;
          return false;
        }
        memcpy(&_nRows, (void *)&(buf[bufPtr]), sizeof(long));

        bufPtr -= sizeof(long);
        if(bufPtr < 0) {
          cerr << __FILE__ << " function:" << __FUNCTION__ \
              << " Error: the buffer has insufficient memory to extract the number of columns" \
              << endl;
          return false;
        }
        memcpy(&_nCols, (void *)&(buf[bufPtr]), sizeof(long));

        _size = _nRows * _nCols;

        buf.erase(buf.begin()+bufPtr, buf.end());

        return true;
      }


    private:
      long _nRows;
      long _nCols;
      long _size;
  };

  /***********************************************
  * CellCoord class                              *
  * Defining the coordinate of a cell:           *
  * Row index and Column index                   *
  ***********************************************/
  class CellCoord {
    public:
      CellCoord()
        :_iRow(pRPL::ERROR_COORD),
         _iCol(pRPL::ERROR_COORD) {}
      CellCoord(long iRow, long iCol)
        :_iRow(iRow),
         _iCol(iCol) {}
      CellCoord(const pRPL::CellCoord &rhs)
        :_iRow(rhs._iRow),
         _iCol(rhs._iCol) {}
      CellCoord(long idx, const pRPL::SpaceDims &dims) {
        fromIdx(idx, dims);
      }
      ~CellCoord() {}

      pRPL::CellCoord& operator=(const pRPL::CellCoord &rhs) {
        if(this != &rhs) {
          _iRow = rhs._iRow;
          _iCol = rhs._iCol;
        }
        return *this;
      }
      void operator()(const pRPL::CellCoord &rhs) {
        _iRow = rhs._iRow;
        _iCol = rhs._iCol;
      }
      void operator()(long iRow, long iCol) {
        _iRow = iRow;
        _iCol = iCol;
      }
      bool operator==(const pRPL::CellCoord &rhs) const {
        return (_iRow == rhs._iRow &&
                _iCol == rhs._iCol);
      }
      bool operator!=(const pRPL::CellCoord &rhs) const {
        return (_iRow != rhs._iRow ||
                _iCol != rhs._iCol);
      }
      pRPL::CellCoord operator+(const pRPL::CellCoord &rhs) const{
        return CellCoord(_iRow + rhs._iRow,
                         _iCol + rhs._iCol);
      }
      pRPL::CellCoord operator-(const pRPL::CellCoord &rhs) const{
        return CellCoord(_iRow - rhs._iRow,
                         _iCol - rhs._iCol);
      }
	  
      void iRow(long idxRow) {
        _iRow = idxRow;
      }
      void iCol(long idxCol) {
        _iCol = idxCol;
      }
      long iRow() const {
        return _iRow;
      }
      long iCol() const {
        return _iCol;
      }

      bool valid() const {
        return(_iRow != pRPL::ERROR_COORD &&
               _iCol != pRPL::ERROR_COORD);
      }
      bool valid(const pRPL::SpaceDims& dims) const {
        return (_iRow >= 0 && _iRow < dims.nRows() &&
                _iCol >= 0 && _iCol < dims.nCols());
      }

	  /***********************************************
    * Public method to determine the direction in
	  * which I am located relative to rhs 
	  ************************************************/
      pRPL::MeshDir locatedTo(const pRPL::CellCoord &rhs) const {
        pRPL::MeshDir dir = pRPL::NONE_DIR;
        if(_iRow < rhs._iRow) {
          if(_iCol < rhs._iCol) dir = pRPL::NORTHWEST_DIR;
          else if(_iCol == rhs._iCol) dir = pRPL::NORTH_DIR;
          else if(_iCol > rhs._iCol) dir = pRPL::NORTHEAST_DIR;
        }
        else if(_iRow == rhs._iRow) {
          if(_iCol < rhs._iCol) dir = pRPL::WEST_DIR;
          else if(_iCol == rhs._iCol) dir = pRPL::NONE_DIR;
          else if(_iCol > rhs._iCol) dir = pRPL::EAST_DIR;
        }
        else if(_iRow > rhs._iRow) {
          if(_iCol < rhs._iCol) dir = pRPL::SOUTHWEST_DIR;
          else if(_iCol == rhs._iCol) dir = pRPL::SOUTH_DIR;
          else if(_iCol > rhs._iCol) dir = pRPL::SOUTHEAST_DIR;
        }
        return dir;
      }
	  
      /***********************************************
      * Public method to return the array index of a
      * given row-column coordinate (i.e., dims).
      * The array index (i.e., location) of an item in a
      * matrix is row_index * num_columns + col_index
      ************************************************/
      long toIdx(const pRPL::SpaceDims& dims) const {
        long idx;
        if(!valid(dims)) {
          idx = pRPL::ERROR_COORD;
        }
        else {
          idx = _iRow * dims.nCols() + _iCol;
        }
        return idx;
      }

      void fromIdx(long idx, const pRPL::SpaceDims& dims) {
        if(!dims.validIdx(idx)) {
          _iRow = pRPL::ERROR_COORD;
          _iCol = pRPL::ERROR_COORD;
        }
        else {
          _iRow = idx / dims.nCols();
          _iCol = idx - _iRow * dims.nCols();
          if(!valid(dims)) {
            _iRow = pRPL::ERROR_COORD;
            _iCol = pRPL::ERROR_COORD;
          }
        }
      }

      void add2Buf(vector<char> &buf) {
        size_t bufPtr = buf.size();
        buf.resize(bufPtr + 2*sizeof(long));

        memcpy((void *)&(buf[bufPtr]), &_iCol, sizeof(long));
        bufPtr += sizeof(long);
        memcpy((void *)&(buf[bufPtr]), &_iRow, sizeof(long));
      }

      bool initFromBuf(vector<char> &buf) {
        _iRow = pRPL::ERROR_COORD;
        _iCol = pRPL::ERROR_COORD;

        int bufPtr = buf.size() - sizeof(long);
        if(bufPtr < 0) {
          cerr << __FILE__ << " function:" << __FUNCTION__ \
               << " Error: the buffer has insufficient memory to extract the row coordinate" \
               << endl;
          return false;
        }
        memcpy(&_iRow, (void *)&(buf[bufPtr]), sizeof(long));

        bufPtr -= sizeof(long);
        if(bufPtr < 0) {
          cerr << __FILE__ << " function:" << __FUNCTION__ \
               << " Error: the buffer has insufficient memory to extract the column coordinate" \
               << endl;
          return false;
        }
        memcpy(&_iCol, (void *)&(buf[bufPtr]), sizeof(long));

        buf.erase(buf.begin()+bufPtr, buf.end());

        return true;
      }

    protected:
      long _iRow;
      long _iCol;
  };

  /***********************************************
  * WeightedCellCoord class                      *
  * Defining the weighted coordinate of a cell:  *
  * Row index, Column index, Weight              *
  ***********************************************/
  class WeightedCellCoord: public pRPL::CellCoord {
    public:
      WeightedCellCoord()
        :pRPL::CellCoord(),
         _weight(1.0) {}
      WeightedCellCoord(long iRow, long iCol,
                        double weight = 1.0)
        :pRPL::CellCoord(iRow, iCol),
         _weight(weight) {}
      WeightedCellCoord(long idx, const pRPL::SpaceDims &dims,
                        double weight = 1.0)
        :pRPL::CellCoord(idx, dims),
         _weight(weight) {}

      WeightedCellCoord(const pRPL::WeightedCellCoord &rhs) {
        pRPL::CellCoord::_iRow = rhs.CellCoord::_iRow;
        pRPL::CellCoord::_iCol = rhs.CellCoord::_iCol;
        _weight = rhs._weight;
      }
      ~WeightedCellCoord() {}
     
      pRPL::WeightedCellCoord& operator=(const pRPL::WeightedCellCoord &rhs) {
        if(this != &rhs) {
          pRPL::CellCoord::_iRow = rhs.CellCoord::_iRow;
          pRPL::CellCoord::_iCol = rhs.CellCoord::_iCol;
          _weight = rhs._weight;
        }
        return *this;
      }
      void operator()(const pRPL::WeightedCellCoord &rhs) {
        pRPL::CellCoord::_iRow = rhs.CellCoord::_iRow;
        pRPL::CellCoord::_iCol = rhs.CellCoord::_iCol;
        _weight = rhs._weight;
      }
      void operator()(long iRow, long iCol,
                      double weight = 1.0) {
        pRPL::CellCoord::_iRow = iRow;
        pRPL::CellCoord::_iCol = iCol;
        _weight = weight;
      }
      bool operator==(const pRPL::WeightedCellCoord &rhs) const {
        return (_iRow == rhs.CellCoord::_iRow &&
                _iCol == rhs.CellCoord::_iCol &&
                std::fabs(float(_weight - rhs._weight)) < EPSINON);
      }
      bool operator!=(const pRPL::WeightedCellCoord &rhs) const {
        return (_iRow != rhs.CellCoord::_iRow ||
                _iCol != rhs.CellCoord::_iCol ||
                std::fabs(float(_weight - rhs._weight)) >= EPSINON);
      }
      void weight(double wVal) {
        _weight = wVal;
      }
      double weight() const {
        return _weight;
      }

      void add2Buf(vector<char> &buf) {
        size_t bufPtr = buf.size();
        buf.resize(bufPtr + sizeof(double));
        memcpy((void *)&(buf[bufPtr]), &_weight, sizeof(double));

        pRPL::CellCoord::add2Buf(buf);
      }

      bool initFromBuf(vector<char> &buf) {
        if(!pRPL::CellCoord::initFromBuf(buf)) {
          return false;
        }

        _weight = 1.0;
        int bufPtr = buf.size() - sizeof(double);
        if(bufPtr < 0) {
          cerr << __FILE__ << " function:" << __FUNCTION__ \
               << " Error: the buffer has insufficient memory to extract the weight" \
               << endl;
          return false;
        }
        memcpy(&_weight, (void *)&(buf[bufPtr]), sizeof(double));

        buf.erase(buf.begin()+bufPtr, buf.end());

        return true;
      }

    protected:
      double _weight;
  };
  
  /***********************************************
  * CoordBR class                                *
  * Defining a bounding rectangle:               *
  * Coordinate of the Northwest Corner Cell;     *
  * Coordinate of the Southeast Corner Cell      *
  ***********************************************/
  class CoordBR {
    public:
      CoordBR() {}
      CoordBR(const pRPL::CellCoord& nw,
              const pRPL::CellCoord& se)
        :_nwCorner(nw),
         _seCorner(se) {}
      CoordBR(const CoordBR& rhs)
        :_nwCorner(rhs._nwCorner),
         _seCorner(rhs._seCorner) {}
      ~CoordBR() {}

      CoordBR& operator=(const pRPL::CoordBR &rhs) {
        if(this != &rhs) {
          _nwCorner = rhs._nwCorner;
          _seCorner = rhs._seCorner;
        }
        return *this;
      }
      bool operator==(const pRPL::CoordBR &rhs) const {
        return (_nwCorner == rhs._nwCorner &&
                _seCorner == rhs._seCorner);
      }
      bool operator!=(const pRPL::CoordBR &rhs) const {
        return (_nwCorner != rhs._nwCorner ||
                _seCorner != rhs._seCorner);
      }
      void nwCorner(const pRPL::CellCoord &nw) {
        _nwCorner = nw;
      }
      void nwCorner(long iRow, long iCol) {
        _nwCorner(iRow, iCol);
      }
      void seCorner(const pRPL::CellCoord &se) {
        _seCorner = se;
      }
      void seCorner(long iRow, long iCol) {
        _seCorner(iRow, iCol);
      }
      
      const pRPL::CellCoord& nwCorner() const {
        return _nwCorner;  
      }
      const pRPL::CellCoord& seCorner() const {
        return _seCorner;  
      }
      
      long minIRow() const {
        return _nwCorner.iRow();
      }
      long minICol() const {
        return _nwCorner.iCol();
      }
      long maxIRow() const {
        return _seCorner.iRow();
      }
      long maxICol() const {
        return _seCorner.iCol();
      }

      bool valid() const {
        return (_nwCorner.valid() && _seCorner.valid() &&
                minIRow() <= maxIRow() && minICol() <= maxICol());
      }
      bool valid(const SpaceDims& dims) const {
        return (_nwCorner.valid(dims) && _seCorner.valid(dims) &&
                valid());
      }
      
      long nRows() const {
        return (valid() ? (maxIRow() - minIRow() + 1):0);
      }
      long nCols() const {
        return (valid() ? (maxICol() - minICol() + 1):0);
      }
      long size() const {
        return (nRows() * nCols());
      }

      bool ifContain(const pRPL::CellCoord &coord) const {
        return (valid() && coord.valid() &&
                coord.iRow() >= minIRow() && coord.iRow() <= maxIRow() &&
                coord.iCol() >= minICol() && coord.iCol() <= maxICol());
      }
      bool ifContain(long idx, const SpaceDims& dims) const {
        return(valid(dims) && ifContain(CellCoord(idx, dims)));
      }
      bool ifContain(const pRPL::CoordBR &rhs) const {
        return (rhs.valid() &&
                ifContain(rhs._nwCorner) &&
                ifContain(rhs._seCorner));
      }
      bool ifWithin(const pRPL::CoordBR &rhs) const {
        return (rhs.ifContain(*this));
      }

      pRPL::CoordBR overlap(const pRPL::CoordBR &rhs) const {
        pRPL::CoordBR olBR;
        if(valid() && rhs.valid()) {
          long iNRow = max(minIRow(), rhs.minIRow());
          long iSRow = min(maxIRow(), rhs.maxIRow());
          long iWCol = max(minICol(), rhs.minICol());
          long iECol = min(maxICol(), rhs.maxICol());
          olBR.nwCorner(iNRow, iWCol);
          olBR.seCorner(iSRow, iECol);
          if(!olBR.valid()) {
            olBR = CoordBR();
          }
        }
        return olBR;
      }

      pRPL::CoordBR merge(const pRPL::CoordBR &rhs) const {
        pRPL::CoordBR mgBR;
        if(valid() && rhs.valid()) {
          long iNRow = min(minIRow(), rhs.minIRow());
          long iSRow = max(maxIRow(), rhs.maxIRow());
          long iWCol = min(minICol(), rhs.minICol());
          long iECol = max(maxICol(), rhs.maxICol());
          mgBR.nwCorner(iNRow, iWCol);
          mgBR.seCorner(iSRow, iECol);
          if(!mgBR.valid()) {
            mgBR = CoordBR();
          }
        }
        return mgBR;
      }
      
      void add2Buf(vector<char> &buf) {
        _seCorner.add2Buf(buf);
        _nwCorner.add2Buf(buf);
      }

      bool initFromBuf(vector<char> &buf) {
        return (_nwCorner.initFromBuf(buf) && _seCorner.initFromBuf(buf));
      }

    private:
      pRPL::CellCoord _nwCorner;
      pRPL::CellCoord _seCorner;
  };

  class GeoCoord {
    public:
    	GeoCoord()
        :_x(pRPL::ERROR_COORD),
         _y(pRPL::ERROR_COORD) {}
      GeoCoord(double x, double y)
        :_x(x),
         _y(y) {}
      GeoCoord(const pRPL::GeoCoord &rhs)
        :_x(rhs._x),
         _y(rhs._y) {}
      ~GeoCoord() {}

      pRPL::GeoCoord& operator=(const pRPL::GeoCoord &rhs) {
        if(this != &rhs) {
          _x = rhs._x;
          _y = rhs._y;
        }
        return *this;
      }
      void operator()(const pRPL::GeoCoord &rhs) {
        _x = rhs._x;
        _y = rhs._y;
      }
      void operator()(double x, double y) {
        _x = x;
        _y = y;
      }
      bool operator==(const pRPL::GeoCoord &rhs) const {
        return (std::fabs(float(_x - rhs._x)) < EPSINON &&
                std::fabs(float(_y - rhs._y)) < EPSINON);
      }
      bool operator!=(const pRPL::GeoCoord &rhs) const {
        return (std::fabs(float(_x - rhs._x)) >=  EPSINON ||
                std::fabs(float(_y - rhs._y)) >= EPSINON);
      }
      pRPL::GeoCoord operator+(const pRPL::GeoCoord &rhs) const{
        return GeoCoord(_x + rhs._x,
                        _y + rhs._y);
      }
      pRPL::GeoCoord operator-(const pRPL::GeoCoord &rhs) const{
        return GeoCoord(_x - rhs._x,
                        _y - rhs._y);
      }
      void x(double xCoord) {
        _x = xCoord;
      }
      void y(double yCoord) {
        _y = yCoord;
      }
      double x() const {
        return _x;
      }
      double y() const {
        return _y;
      }
      
      void add2Buf(vector<char> &buf) {
        size_t bufPtr = buf.size();
        buf.resize(bufPtr + 2*sizeof(double));

        memcpy((void *)&(buf[bufPtr]), &_y, sizeof(double));
        bufPtr += sizeof(double);
        memcpy((void *)&(buf[bufPtr]), &_x, sizeof(double));
      }

      bool initFromBuf(vector<char> &buf)  {
        _x = pRPL::ERROR_COORD;
        _y = pRPL::ERROR_COORD;

        int bufPtr = buf.size() - sizeof(double);
        if(bufPtr < 0) {
          cerr << __FILE__ << " function:" << __FUNCTION__ \
              << " Error: the buffer has insufficient memory to extract the X coordinate" \
              << endl;
          return false;
        }
        memcpy(&_x, (void *)&(buf[bufPtr]), sizeof(double));

        bufPtr -= sizeof(double);
        if(bufPtr < 0) {
          cerr << __FILE__ << " function:" << __FUNCTION__ \
              << " Error: the buffer has insufficient memory to extract the Y coordinate" \
              << endl;
          return false;
        }
        memcpy(&_y, (void *)&(buf[bufPtr]), sizeof(double));

        buf.erase(buf.begin()+bufPtr, buf.end());

        return true;
      }

    private:
      double _x;
      double _y;
  };
  
  /***********************************************
  * IO-Stream Operators for SpaceDims,           *
  * CellCoord, CoordBR, and GeoCoord             *
  ***********************************************/
  ostream& operator<<(ostream &os, const pRPL::SpaceDims &dims);
  istream& operator>>(istream &is, pRPL::SpaceDims &dims);
  ostream& operator<<(ostream &os, const pRPL::CellCoord &coord);
  istream& operator>>(istream &is, pRPL::CellCoord &coord);
  ostream& operator<<(ostream &os, const pRPL::WeightedCellCoord &coord);
  istream& operator>>(istream &is, pRPL::WeightedCellCoord &coord);
  ostream& operator<<(ostream &os, const pRPL::CoordBR &br);
  ostream& operator<<(ostream &os, const pRPL::GeoCoord &gcoord);
  istream& operator>>(istream &is, pRPL::GeoCoord &gcoord);
};


#endif
