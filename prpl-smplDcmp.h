#ifndef PRPL_SMPLDCMP_H
#define PRPL_SMPLDCMP_H

#include "prpl-basicTypes.h"
#include "prpl-neighborhood.h"
#include "prpl-subCellspaceInfo.h"

namespace pRPL {
  class SmplDcmp{
    public:
      SmplDcmp() {}
      ~SmplDcmp() {}

      bool rowDcmp(list<pRPL::SubCellspaceInfo> &lSubspcInfos,
                   const pRPL::SpaceDims &glbDims,
                   int nSubspcs,
                   const pRPL::Neighborhood *pNbrhd = NULL) const;
      bool colDcmp(list<pRPL::SubCellspaceInfo> &lSubspcInfos,
                   const pRPL::SpaceDims &glbDims,
                   int nSubspcs,
                   const pRPL::Neighborhood *pNbrhd = NULL) const;
      bool blkDcmp(list<pRPL::SubCellspaceInfo> &lSubspcInfos,
                   const pRPL::SpaceDims &glbDims,
                   int nRowSubspcs,
                   int nColSubspcs = 1,
                   const pRPL::Neighborhood *pNbrhd = NULL) const;

    private:
      bool _checkDcmpPrmtrs(list<pRPL::SubCellspaceInfo> &lSubspcInfos,
                            const pRPL::SpaceDims &glbDims,
                            int nRowSubspcs, int nColSubspcs,
                            const pRPL::Neighborhood *pNbrhd) const;

      void _smpl1DDcmp(long &subBegin, long &subEnd,
                       long glbBegin, long glbEnd,
                       int nSubspcs, int iSubspc) const;
  };
};


#endif /* SMPLDCMP_H */
