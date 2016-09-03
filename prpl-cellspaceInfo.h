#ifndef PRPL_CELLSPACEINFO_H
#define PRPL_CELLSPACEINFO_H

#include "prpl-basicTypes.h"
#include "prpl-cellspaceGeoinfo.h"
#include "prpl-neighborhood.h"

namespace pRPL {
  class CellspaceInfo {
    public:
      CellspaceInfo();
      CellspaceInfo(const pRPL::CellspaceInfo &rhs);
      CellspaceInfo(const pRPL::SpaceDims &dims,
                    const char *aTypeName,
                    size_t typeSize,
                    const pRPL::CellspaceGeoinfo *pGeoinfo = NULL,
                    long tileSize = TILEWIDTH);
      ~CellspaceInfo();

      pRPL::CellspaceInfo& operator=(const pRPL::CellspaceInfo &rhs);
      bool operator==(const pRPL::CellspaceInfo &rhs);
      bool operator!=(const pRPL::CellspaceInfo &rhs);

      bool init(const pRPL::SpaceDims &dims,
                const char *aTypeName,
                size_t typeSize,
                const pRPL::CellspaceGeoinfo *pGeoinfo = NULL,
                long tileSize = TILEWIDTH);

      bool initByGDAL(GDALDataset *pDataset,
                      int iBand,
                      bool warning = true);
     bool initByPGTIOL(PGTIOLDataset *pDataset,
                       int iBand,
                       bool warning = true);

      void clear();

      void add2Buf(vector<char> &buf);
      bool initFromBuf(vector<char> &buf);

      bool isEmpty(bool warning = false) const;
      void dims(const pRPL::SpaceDims& cellDims);
      const pRPL::SpaceDims& dims() const;
      long nRows() const;
      long nCols() const;

      long size() const;
      const char* dataType() const;
      size_t dataSize() const;

      void georeference(const pRPL::CellspaceGeoinfo &geoinfo);
      const pRPL::CellspaceGeoinfo* georeference() const;
      bool isGeoreferenced(bool warning = false) const;

      template<class CmprType>
      bool isDataType(bool warning = false) const;

      template<class RtrType>
      RtrType& getNoDataVal(bool warning = true) const;
      template<class RtrType>
      const RtrType getNoDataValAs(bool warning = true) const;

      template<class ElemType>
      bool setNoDataVal(const ElemType &val,
                        bool warning = true);
      template<class ElemType>
      bool setNoDataValAs(const ElemType &val,
                          bool warning = true);

      long tileSize() const;
      void setTileSize(long size);

      /* Coordinate */
      long coord2idx(const pRPL::CellCoord &coord) const;
      long coord2idx(long iRow, long iCol) const;
      pRPL::CellCoord idx2coord(long idx) const;
      bool validCoord(const pRPL::CellCoord &coord,
                      bool warning = true) const;
      bool validCoord(long iRow, long iCol,
                      bool warning = true) const;
      bool validIdx(long idx,
                    bool warning = true) const;
      bool validBR(const pRPL::CoordBR &rectangle,
                   bool warning = true) const;

      pRPL::GeoCoord cellCoord2geoCoord(const pRPL::CellCoord &cCoord) const;
      pRPL::GeoCoord cellCoord2geoCoord(long iRow, long iCol) const;
      pRPL::CellCoord geoCoord2cellCoord(const pRPL::GeoCoord &gCoord) const;
      pRPL::CellCoord geoCoord2cellCoord(double x, double y) const;

      /* Neighborhood operations */
      bool calcWorkBR(pRPL::CoordBR *pWorkBR,
                      const pRPL::Neighborhood *pNbrhd,
                      bool warning = true) const;

      bool validNbrhd(const pRPL::Neighborhood &nbrhd,
                      const pRPL::CellCoord &ctrCoord,
                      bool warning = true) const;
      bool validNbrhd(const pRPL::Neighborhood &nbrhd,
                      long ctrRow, long ctrCol,
                      bool warning = true) const;
      bool validNbrhd(const pRPL::Neighborhood &nbrhd,
                      long idx,
                      bool warning = true) const;

    private:
      pRPL::SpaceDims _dims;
      string _dTypeName;
      size_t _dTypeSize;
      void *_pNoDataVal;
      pRPL::CellspaceGeoinfo *_pGeoInfo;
      long _tileSize;
  };
};

template<class CmprType>
bool pRPL::CellspaceInfo::
isDataType(bool warning) const {
  bool valid = true;
  if(_dTypeName.compare(typeid(CmprType).name()) || _dTypeSize != sizeof(CmprType)) {
    if(warning) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: Cellspace element data type (" << _dTypeName \
           << ") does NOT match requested data type (" << typeid(CmprType).name() \
           << ")" << endl;
    }
    valid = false;
  }
  return valid;
}

template<class RtrType>
RtrType& pRPL::CellspaceInfo::
getNoDataVal(bool warning) const {
  RtrType *pNoData = NULL;
  if(isDataType<RtrType>(warning) && _pNoDataVal != NULL) {
    pNoData = (RtrType*)_pNoDataVal;
  }
  return *pNoData;
}

template<class RtrType>
const RtrType pRPL::CellspaceInfo::
getNoDataValAs(bool warning) const {
  RtrType val = (RtrType)pRPL::ERROR_VAL;

  if(_pNoDataVal != NULL) {
    if(isDataType<bool>()) {
      val = (RtrType)(getNoDataVal<bool>());
    }
    else if(isDataType<unsigned char>()) {
      val = (RtrType)(getNoDataVal<unsigned char>());
    }
    else if(isDataType<char>()) {
      val = (RtrType)(getNoDataVal<char>());
    }
    else if(isDataType<unsigned short>()) {
      val = (RtrType)(getNoDataVal<unsigned short>());
    }
    else if(isDataType<short>()) {
      val = (RtrType)(getNoDataVal<short>());
    }
    else if(isDataType<unsigned int>()) {
      val = (RtrType)(getNoDataVal<unsigned int>());
    }
    else if(isDataType<int>()) {
      val = (RtrType)(getNoDataVal<int>());
    }
    else if(isDataType<unsigned long>()) {
      val = (RtrType)(getNoDataVal<unsigned long>());
    }
    else if(isDataType<long>()) {
      val = (RtrType)(getNoDataVal<long>());
    }
    else if(isDataType<float>()) {
      val = (RtrType)(getNoDataVal<float>());
    }
    else if(isDataType<double>()) {
      val = (RtrType)(getNoDataVal<double>());
    }
    else if(isDataType<long double>()) {
      val = (RtrType)(getNoDataVal<long double>());
    }
    else {
      if(warning) {
        cerr << __FILE__ << " function:" << __FUNCTION__ \
            << " Error: unable to convert Cellspace data type (" \
            << dataType() << ") to output data type (" << typeid(RtrType).name() \
            << "). Returned an INVALID value. "<< endl;
      }
    }
  } // end -- if(_pNoDataVal != NULL)
  return val;
}

template<class ElemType>
bool pRPL::CellspaceInfo::
setNoDataVal(const ElemType &val,
             bool warning) {
  bool done = false;
  if(isDataType<ElemType>(warning) && _pNoDataVal != NULL) {
    (ElemType&)(*((ElemType*)_pNoDataVal)) = val;
    done = true;
  }
  return done;
}

template<class ElemType>
bool pRPL::CellspaceInfo::
setNoDataValAs(const ElemType &val,
               bool warning) {
  bool done = false;
  if(!isEmpty(warning) && _pNoDataVal != NULL) {
    if(isDataType<bool>()) {
      done = setNoDataVal(static_cast<bool>(val), warning);
    }
    else if(isDataType<unsigned char>()) {
      done = setNoDataVal(static_cast<unsigned char>(val), warning);
    }
    else if(isDataType<char>()) {
      done = setNoDataVal(static_cast<char>(val), warning);
    }
    else if(isDataType<unsigned short>()) {
      done = setNoDataVal(static_cast<unsigned short>(val), warning);
    }
    else if(isDataType<short>()) {
      done = setNoDataVal(static_cast<short>(val), warning);
    }
    else if(isDataType<unsigned int>()) {
      done = setNoDataVal(static_cast<unsigned int>(val), warning);
    }
    else if(isDataType<int>()) {
      done = setNoDataVal(static_cast<int>(val), warning);
    }
    else if(isDataType<unsigned long>()) {
      done = setNoDataVal(static_cast<unsigned long>(val), warning);
    }
    else if(isDataType<long>()) {
      done = setNoDataVal(static_cast<long>(val), warning);
    }
    else if(isDataType<float>()) {
      done = setNoDataVal(static_cast<float>(val), warning);
    }
    else if(isDataType<double>()) {
      done = setNoDataVal(static_cast<double>(val), warning);
    }
    else if(isDataType<long double>()) {
      done = setNoDataVal(static_cast<long double>(val), warning);
    }
    else {
      if(warning) {
        cerr << __FILE__ << " function:" << __FUNCTION__ \
             << " Error: unable to convert input data type (" \
             << typeid(ElemType).name() << ") to Cellspace data type (" \
             << dataType() << ")" << endl;
        done = false;
      }
    }
  } // end -- if(!isEmpty(warning) && _pNoDataVal != NULL)

  return done;
}
#endif /* CELLSPACEINFO_H */
