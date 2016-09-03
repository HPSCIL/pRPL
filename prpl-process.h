#ifndef PRPL_PROCESS_H
#define PRPL_PROCESS_H
#define MPICH_SKIP_MPICXX
#include "prpl-basicTypes.h"
#include "mpi.h"

namespace pRPL {
  class Process {
    public:
      Process();
      ~Process() {}

      bool initialized() const;
      bool active() const;
      bool init(int argc = 0,
                char *argv[] = NULL);
      void abort() const;
      void finalize() const;
      void sync() const;

      bool set(MPI_Comm &comm,
               bool hasWriter = false,
               int groupID = -1);
      bool grouping(int nGroups,
                    bool incldMaster,
                    Process *pGrpedPrc,
                    Process *pGrpMaster = NULL) const;

      const MPI_Comm& comm() const;
      int id() const;
      const char* processorName() const;
      int groupID() const;
      int nProcesses() const;
      bool hasWriter() const;
      bool isMaster() const;
      bool isWriter() const;
      const pRPL::IntVect allPrcIDs(bool incldMaster = false,
                                    bool incldWriter = false) const;

    private:
      MPI_Comm _comm;
      int _id;
      int _grpID;
      int _nTotalPrcs;
      bool _hasWriter;
      string _prcrName;
  };
};


#endif /* PROCESS_H */
