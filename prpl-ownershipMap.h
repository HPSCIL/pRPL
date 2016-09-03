#ifndef PRPL_OWNERSHIPMAP_H
#define PRPL_OWNERSHIPMAP_H

#include "prpl-basicTypes.h"

namespace pRPL {
  class OwnershipMap {
    public:
      OwnershipMap();
      ~OwnershipMap() {}

      void init(const pRPL::IntVect &vSubspcIDs,
                const pRPL::IntVect &vPrcIDs);
      void clearMapping();

      void mapping(pRPL::MappingMethod mapMethod = pRPL::CYLC_MAP,
                   int nSubspcs2Map = 0);
      int mappingNextTo(int prcID); // return the SubCellspace's ID
      bool mappingTo(int subspcID,
                     int prcID);
      bool allMapped() const;

      const pRPL::IntVect* subspcIDsOnPrc(int prcID) const;
      pRPL::IntVect mappedSubspcIDs() const;
      int findPrcBySubspcID(int subspcID) const;
      bool findSubspcIDOnPrc(int subspcID, int prcID) const;

      bool map2buf(vector<char> &buf) const;
      bool buf2map(const vector<char> &buf);

    private:
      pRPL::IntVect _vSubspcIDs;
      map<int, pRPL::IntVect> _mPrcSubspcs;

      pRPL::IntVect::iterator _itrSubspc2Map;

    friend ostream& operator<<(ostream &os, const pRPL::OwnershipMap &mOwnerships);
  };

  ostream& operator<<(ostream &os, const pRPL::OwnershipMap &mOwnerships);
};


#endif /* OwnershipMap_H */
