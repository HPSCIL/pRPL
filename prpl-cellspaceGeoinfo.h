#ifndef PRPL_CELLSPACEGEOINFO_H
#define PRPL_CELLSPACEGEOINFO_H

#include "prpl-basicTypes.h"
#include "prpl-subCellspaceInfo.h"

namespace pRPL {
  class CellspaceGeoinfo {
    public:
      CellspaceGeoinfo() {}
      CellspaceGeoinfo(const pRPL::GeoCoord &nwCorner,
                       const pRPL::GeoCoord &cellSize,
                       const string &projection);
      CellspaceGeoinfo(const pRPL::CellspaceGeoinfo &rhs);
      ~CellspaceGeoinfo() {}

      pRPL::CellspaceGeoinfo& operator=(const pRPL::CellspaceGeoinfo &rhs);
      bool operator==(const pRPL::CellspaceGeoinfo &rhs) const;

      void add2Buf(vector<char> &buf);
      bool initFromBuf(vector<char> &buf);

      void nwCorner(const pRPL::GeoCoord &nw);
      void nwCorner(double x, double y);
      void cellSize(const pRPL::GeoCoord &size);
      void cellSize(double xSize, double ySize);
      void projection(const string &proj);

      const pRPL::GeoCoord& nwCorner() const;
      const pRPL::GeoCoord& cellSize() const;

      const string& projection() const;
      void geoTransform(double aGeoTransform[6]) const;

      pRPL::GeoCoord cellCoord2geoCoord(const pRPL::CellCoord &cCoord) const;
      pRPL::GeoCoord cellCoord2geoCoord(long iRow, long iCol) const;
      pRPL::CellCoord geoCoord2cellCoord(const pRPL::GeoCoord &gCoord) const;
      pRPL::CellCoord geoCoord2cellCoord(double x, double y) const;

	  /*--------GDAL------------*/
      bool initByGDAL(GDALDataset *pDataset,
                      bool warning = true);
	  /*---------PGTIOL--------*/
	  bool initByPGTIOL(PGTIOLDataset *pDataset,
                      bool warning = true);
	  
      pRPL::CellspaceGeoinfo subSpaceGeoinfo(const pRPL::SubCellspaceInfo &subInfo) const;

    private:
      pRPL::GeoCoord _nwCorner;
      pRPL::GeoCoord _cellSize;
      string _projection;
  };  
};

#endif
