#ifndef PRPL_EXCHANGEMAP_H
#define PRPL_EXCHANGEMAP_H

#include "prpl-basicTypes.h"

namespace pRPL {
  typedef struct {
    int subspcGlbID;
    int iDirection;
    int iNeighbor;
  } ExchangeNode;

  class ExchangeMap: public map<int, vector<pRPL::ExchangeNode> > {
    public:
      bool add(int subspcID, int iDir, int iNbr, int prcID);
  };

};

#endif
