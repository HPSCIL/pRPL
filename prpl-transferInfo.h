#ifndef PRPL_TRANSFERINFO_H
#define PRPL_TRANSFERINFO_H

#include "prpl-basicTypes.h"
#include "prpl-transition.h"

namespace pRPL {
  class TransferInfo {
    public:
      TransferInfo();
      TransferInfo(int fromPrcID,
                   int toPrcID,
                   const string &lyrName,
                   int subspcGlbID = pRPL::ERROR_ID,
                   pRPL::TransferDataType tDataType = pRPL::SPACE_TRANSFDATA);
      ~TransferInfo() {}

      int fromPrcID() const;
      int toPrcID() const;
      const string& lyrName() const;
      int subspcGlbID() const;
      pRPL::TransferDataType transfDataType() const;
      bool completed() const;

      void complete();
      void restart();

    protected:
      int _fromPrcID;
      int _toPrcID;
      string _lyrName;
      int _subspcGlbID;
      pRPL::TransferDataType _tDataType; // SPACE_TRANSF or CELL_TRANSF
      bool _cmplt;
  };

  class TransferInfoVect: public vector<pRPL::TransferInfo> {
    public:
      const pRPL::TransferInfo* findInfo(const string &lyrName,
                                         int subspcGlbID) const;

      /* Master checks if a Layer's Cellspace has been broadcasted */
      bool checkBcastCellspace(int prcID,
                               const string &lyrName) const;

  };
};

#endif

