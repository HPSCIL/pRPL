#ifndef PRPL_CELLSTREAM_H
#define PRPL_CELLSTREAM_H

#include "prpl-basicTypes.h"

namespace pRPL {
  class CellStream {
    public:
      CellStream();
      ~CellStream();

      pRPL::CellStream& operator=(const pRPL::CellStream &rhs);

      void clear();

      size_t size() const;
      int getTotalCellCount() const;
      vector<int>& getCellCounts();
      const vector<int>& getCellCounts() const;
      const vector<pair<unsigned int, unsigned int> >& getLayerInfos() const;

      bool resize();
      void* getStream();
      const void* getStream() const;

      bool addLayer(int lyrID,
                    unsigned int typeSize);
      bool addCell(long cellIdx,
                   const void *aCellVal);

      unsigned int getNumLayers() const;
      unsigned int getNumCellsOnLayer(int lyrID) const;
      unsigned int getTypeSizeOnLayer(int lyrID) const;
      void* getCellsOnLayer(int lyrID);

    private:
      vector<pair<unsigned int, unsigned int> >::const_iterator _getLayerInfo(int lyrID) const;

    private:
      size_t _size;
      void *_aStream; // pair<long, ANY_TYPE>: Cell index, Cell value
      vector<int> _vCellCounts; // numbers of cells on Layers
      vector<pair<unsigned int, unsigned int> > _vLayerInfos; // Layer index, data type size
  };
};

#endif
