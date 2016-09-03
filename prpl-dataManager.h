#ifndef PRPL_DATAMANAGER_H
#define PRPL_DATAMANAGER_H

#include <time.h>
#include "prpl-basicTypes.h"
#include "prpl-cellspace.h"
#include "prpl-layer.h"
#include "prpl-neighborhood.h"
#include "prpl-transition.h"
#include "prpl-smplDcmp.h"
#include "prpl-ownershipMap.h"
#include "prpl-exchangeMap.h"
#include "prpl-process.h"
#include "prpl-transferInfo.h"

namespace pRPL {
  class DataManager {
    public:
      /* Constructor and destructor */
      DataManager();
      ~DataManager();
      
      /* MPI Process */
      bool initMPI(MPI_Comm comm = MPI_COMM_WORLD,
                   bool hasWriter = false);
      void finalizeMPI();
      pRPL::Process& mpiPrc();

      /* Layers */
      int nLayers() const;
      pRPL::Layer* addLayer(const char *aLyrName);

      /*ByGDAL----*/
      pRPL::Layer* addLayerByGDAL(const char *aLyrName,
                                  const char *aGdalFileName,
                                  int iBand,
                                  bool pReading = true);
      bool createLayerGDAL(const char *aLyrName,
                           const char *aGdalFileName,
                           const char *aGdalFormat,
                           char **aGdalOptions = NULL);
      void closeGDAL();

      /*ByPGTIOL--hsj-*/
      pRPL::Layer* addLayerByPGTIOL(const char *aLyrName,
                                    const char *aPgtiolFileName,
                                    int iBand,
                                    bool pReading = true);
      bool createLayerPGTIOL(const char *aLyrName,
                             const char *aPgtiolFileName,
                             char **aPgtiolOptions = NULL);
      void closePGTIOL();
      void closeDatasets();

      bool rmvLayerByIdx(int lyrID,
                         bool warning = true);
      bool rmvLayerByName(const char *aLyrName,
                          bool warning = true);
      void clearLayers();
      pRPL::Layer* getLayerByIdx(int lyrID,
                                 bool wrning = true);
      const pRPL::Layer* getLayerByIdx(int lyrID,
                                       bool wrning = true) const ;
      int getLayerIdxByName(const char *aLyrName,
                            bool warning = true) const ;
      pRPL::Layer* getLayerByName(const char *aLyrName,
                                  bool warning = true);
      const pRPL::Layer* getLayerByName(const char *aLyrName,
                                        bool warning = true) const;
      const char* getLayerNameByIdx(int lyrID,
                                    bool warning = true) const;
      bool beginReadingLayer(const char *aLyrName,
                             pRPL::ReadingOption readOpt);
      void finishReadingLayers(pRPL::ReadingOption readOpt);

      /* Neighborhoods */
      int nNbrhds() const;
      pRPL::Neighborhood* addNbrhd(const char *aNbrhdName);
      pRPL::Neighborhood* addNbrhd(const pRPL::Neighborhood &rhs);
      bool rmvNbrhdByIdx(int nbrhdID,
                         bool warning = true);
      bool rmvNbrhdByName(const char *aNbrhdName,
                          bool warning = true);
      void clearNbrhds();
      pRPL::Neighborhood* getNbrhdByIdx(int nbrhdID = 0,
                                        bool warning = true);
      const pRPL::Neighborhood* getNbrhdByIdx(int nbrhdID = 0,
                                              bool warning = true) const;
      int getNbrhdIdxByName(const char *aNbrhdName,
                            bool warning = true);
      pRPL::Neighborhood* getNbrhdByName(const char *aNbrhdName,
                                         bool warning = true);
      const pRPL::Neighborhood* getNbrhdByName(const char *aNbrhdName,
                                               bool warning = true) const;
      
      /* Decomposition */
      bool dcmpLayer(int lyrID,
                     int nbrhdID,
                     int nRowSubspcs,
                     int nColSubspcs);
      bool dcmpLayer(const char *aLyrName,
                     const char *aNbrhdName,
                     int nRowSubspcs,
                     int nColSubspcs);

      bool propagateDcmp(int fromLyrID,
                         int toLyrID);
      bool propagateDcmp(int fromLyrID,
                         const pRPL::IntVect &vToLyrIDs);
      bool propagateDcmp(const string &fromLyrName,
                         const string &toLyrName);
      bool propagateDcmp(const string &fromLyrName,
                         const vector<string> &vToLyrNames);

      bool dcmpLayers(const pRPL::IntVect &vLyrIDs,
                      int nbrhdID,
                      int nRowSubspcs,
                      int nColSubspcs);
      bool dcmpLayers(const vector<string> &vLyrNames,
                      const char *aNbrhdName,
                      int nRowSubspcs,
                      int nColSubspcs);
      bool dcmpLayers(pRPL::Transition &trans,
                      int nRowSubspcs,
                      int nColSubspcs);
      bool dcmpAllLayers(const char *aNbrhdName,
                         int nRowSubspcs,
                         int nColSubspcs);

      const pRPL::OwnershipMap& ownershipMap() const;
      bool syncOwnershipMap();

      /* Evaluate - Dynamic Tasking (Task Farm) */
      bool initTaskFarm(pRPL::Transition &trans,
                        pRPL::MappingMethod mapMethod = pRPL::CYLC_MAP,
                        int nSubspcs2Map = 0,
                        pRPL::ReadingOption readOpt = pRPL::PARA_READING);
      pRPL::EvaluateReturn evaluate_TF(pRPL::EvaluateType evalType,
                                       pRPL::Transition &trans,
                                       pRPL::ReadingOption readOpt,
                                       pRPL::WritingOption writeOpt,
                                       bool ifInitNoData = false,
                                       bool ifSyncOwnershipMap = false,
                                       const pRPL::LongVect *pvGlbIdxs = NULL);

      /* Evaluate - Static Tasking */
      bool initStaticTask(pRPL::Transition &trans,
                          pRPL::MappingMethod mapMethod = pRPL::CYLC_MAP,
                          pRPL::ReadingOption readOpt = pRPL::PARA_READING);
      pRPL::EvaluateReturn evaluate_ST(pRPL::EvaluateType evalType,
                                       pRPL::Transition &trans,
                                       pRPL::WritingOption writeOpt,
                                       bool ifInitNoData = false,
                                       const pRPL::LongVect *pvGlbIdxs = NULL);

      /* Utilities */
      bool mergeTmpSubspcByGDAL(pRPL::Transition &trans,
                                int subspcGlbID);
      bool mergeAllTmpSubspcsByGDAL(pRPL::Transition &trans);

    protected:
      bool _initReadInData(pRPL::Transition &trans,
                           pRPL::ReadingOption readOpt = pRPL::PARA_READING);

      bool _toMsgTag(int &msgTag,
                     int lyrID,
                     int subspcID,
                     int maxSubspcID);

      bool _iBcastCellspaceBegin(const string &lyrName);
      bool _iTransfSubspaceBegin(int fromPrcID,
                                 int toPrcID,
                                 const string &lyrName,
                                 int subspcGlbID);

      // returns the number of completed non-blocking (Sub)Cellspace communications
      // returns MPI_UNDEFINED if there is no active request
      int _iTransfSpaceTest(pRPL::IntVect &vCmpltIDs,
                            pRPL::TransferType transfType,
                            bool delSentSpaces = false,
                            bool writeRecvSpaces = false,
                            bool delRecvSpaces = false);

      void _clearTransfSpaceReqs();

      int _checkTransInDataReady(pRPL::Transition &trans,
                                 pRPL::ReadingOption readOpt);
      bool _initTransOutData(pRPL::Transition &trans,
                             int subspcGlbID,
                             bool ifInitNoData = false);

      bool _setTransCellspaces(pRPL::Transition &trans,
                               int subspcGlbID = pRPL::ERROR_ID);
      bool _setTransNbrhd(pRPL::Transition &trans);

      bool _initEvaluate(pRPL::Transition &trans,
                         int subspcGlbID,
                         pRPL::EvaluateBR br2Eval = pRPL::EVAL_WORKBR,
                         bool ifInitNoData = false);
      bool _fnlzEvaluate(pRPL::Transition &trans,
                         int subspcGlbID);

      bool _workerReadSubspcs(pRPL::Transition &trans,
                              int subspcGlbID,
                              pRPL::ReadingOption readOpt);
      bool _workerWriteSubspcs(pRPL::Transition &trans,
                               int subspcGlbID,
                               pRPL::WritingOption writeOpt);

      bool _dymWriteSubspcs(pRPL::Transition &trans,
                            int subspcGlbID,
                            pRPL::WritingOption writeOpt);
      bool _finalWriteSubspcs(pRPL::Transition &trans,
                              pRPL::WritingOption writeOpt);

      pRPL::EvaluateReturn _evalNone(pRPL::Transition &trans,
                                     int subspcGlbID = pRPL::ERROR_ID,
                                     bool ifInitNoData = false);
      pRPL::EvaluateReturn _evalAll(pRPL::Transition &trans,
                                    int subspcGlbID = pRPL::ERROR_ID,
                                    pRPL::EvaluateBR br2Eval = pRPL::EVAL_WORKBR,
                                    bool ifInitNoData = false);
      pRPL::EvaluateReturn _evalRandomly(pRPL::Transition &trans,
                                         int subspcGlbID = pRPL::ERROR_ID,
                                         bool ifInitNoData = false);
      pRPL::EvaluateReturn _evalSelected(pRPL::Transition &trans,
                                         const pRPL::LongVect &vGlbIdxs,
                                         int subspcGlbID = pRPL::ERROR_ID,
                                         pRPL::EvaluateBR br2Eval = pRPL::EVAL_WORKBR,
                                         bool ifInitNoData = false);

      pRPL::EvaluateReturn _masterTF(pRPL::Transition &trans,
                                     pRPL::ReadingOption readOpt,
                                     pRPL::WritingOption writeOpt);
      pRPL::EvaluateReturn _writerTF(pRPL::Transition &trans,
                                     pRPL::WritingOption writeOpt);
      pRPL::EvaluateReturn _workerTF(pRPL::EvaluateType evalType,
                                     pRPL::Transition &trans,
                                     pRPL::ReadingOption readOpt,
                                     pRPL::WritingOption writeOpt,
                                     bool ifInitNoData = false,
                                     const pRPL::LongVect *pvGlbIdxs = NULL);

      pRPL::EvaluateReturn _masterST(pRPL::Transition &trans,
                                     pRPL::WritingOption writeOpt);
      pRPL::EvaluateReturn _writerST(pRPL::Transition &trans,
                                     pRPL::WritingOption writeOpt);
      pRPL::EvaluateReturn _workerST(pRPL::EvaluateType evalType,
                                     pRPL::Transition &trans,
                                     pRPL::WritingOption writeOpt,
                                     bool ifInitNoData = false,
                                     const pRPL::LongVect *pvGlbIdxs = NULL);

      bool _calcExchangeBRs(pRPL::Transition &trans);
      bool _calcExchangeRoutes(pRPL::Transition &trans);

      bool _loadCellStream(pRPL::CellStream &cells2load);
      bool _makeCellStream(pRPL::Transition &trans);
      bool _iExchangeBegin();
      bool _iExchangeEnd();

    private:
      pRPL::Process _prc;
      int _nActiveWorkers;

      list<pRPL::Layer> _lLayers;
      list<pRPL::Neighborhood> _lNbrhds;
      pRPL::OwnershipMap _mOwnerships;

      vector<MPI_Request> _vBcstSpcReqs;
      pRPL::TransferInfoVect _vBcstSpcInfos;
      vector<MPI_Request> _vSendSpcReqs;
      pRPL::TransferInfoVect _vSendSpcInfos;
      vector<MPI_Request> _vRecvSpcReqs;
      pRPL::TransferInfoVect _vRecvSpcInfos;
      pRPL::IntVect _vEvaledSubspcIDs;

      pRPL::ExchangeMap _mSendRoutes;
      pRPL::ExchangeMap _mRecvRoutes;
      map<int, pRPL::CellStream> _mSendCells;
      map<int, pRPL::CellStream> _mRecvCells;
      vector<MPI_Request> _vSendCellReqs;
      vector<MPI_Request> _vRecvCellReqs;
  };
};

#endif
