#include "prpl-exchangeMap.h"

bool pRPL::ExchangeMap::
add(int subspcID, int iDir, int iNbr, int prcID) {
  if(prcID < 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid process ID (" << prcID \
         << ")" << endl;
    return false;
  }
  if(subspcID < 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid SubCellspace ID (" << subspcID \
         << ")" << endl;
    return false;
  }
  if(iDir < 0 || iNbr < 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid direction (" << iDir << ")," \
         << " or invalid neighbor ID (" << iNbr << ")" \
         << endl;
    return false;
  }

  pRPL::ExchangeNode newNode = {subspcID, iDir, iNbr};
  operator[](prcID).push_back(newNode);
  //cout << subspcID << ", " << iDir << ", " << iNbr << " ---- " << prcID << endl;
  return true;
}
