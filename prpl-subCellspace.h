#ifndef PRPL_SUBCELLSPACE_H
#define PRPL_SUBCELLSPACE_H

#include "prpl-cellspace.h"
#include "prpl-subCellspaceInfo.h"
#include "prpl-cellStream.h"

namespace pRPL {
  class SubCellspace: public pRPL::Cellspace {
    public:
      /* Constructors and Destructor */
      SubCellspace();
      SubCellspace(pRPL::CellspaceInfo *pInfo,
                   pRPL::SubCellspaceInfo *pSubInfo);
      SubCellspace(const pRPL::SubCellspace &rhs);
      ~SubCellspace();

      /* Operators */
      pRPL::SubCellspace& operator=(const pRPL::SubCellspace &rhs);

      pRPL::SubCellspaceInfo* subInfo();
      const pRPL::SubCellspaceInfo* subInfo() const;

      bool add2UpdtStream(pRPL::CellStream &updtStream,
                          int iDir,
                          int iNbr = 0) const;
      bool loadUpdtStream(const pRPL::CellStream &updtStream);

    protected:
      pRPL::SubCellspaceInfo *_pSubInfo;
  };
};

#endif
