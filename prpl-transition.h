#ifndef PRPL_TRANSITION_H
#define PRPL_TRANSITION_H

#include "prpl-basicTypes.h"
#include "prpl-cellspace.h"
#include "prpl-neighborhood.h"
#include "prpl-subCellspace.h"

namespace pRPL {
  class Transition {
    public:
      /* Constructor and destructor */
      Transition(bool onlyUpdtCtrCell = true,
                 bool needExchange = false,
                 bool edgesFirst = false);
      virtual ~Transition() {}

      /* Layers */
      void addInputLyr(const char *aInLyrName,
                       bool isPrimeLyr = false);
      void addOutputLyr(const char *aOutLyrName,
                        bool isPrimeLyr = true);
      bool setLyrsByNames(vector<string> *pvInLyrNames,
                          vector<string> *pvOutLyrNames,
                          string &primeLyrName);
      const vector<string>& getInLyrNames() const;
      const vector<string>& getOutLyrNames() const;
      const string& getPrimeLyrName() const;
      bool isInLyr(const string &lyrName) const;
      bool isOutLyr(const string &lyrName) const;
      bool isPrimeLyr(const string &lyrName) const;
      void clearLyrSettings();

      bool setCellspace(const string &lyrName,
                        pRPL::Cellspace *pCellspc);
      void clearCellspaces();

      pRPL::Cellspace* getCellspaceByLyrName(const string &lyrName);
      const pRPL::Cellspace* getCellspaceByLyrName(const string &lyrName) const;

      void setUpdateTracking(bool toTrack);
      void clearUpdateTracks();

      /* Neighborhood */
      void setNbrhdByName(const char *aNbrhdName);
      const string& getNbrhdName() const;
      void clearNbrhdSetting();

      void setNbrhd(Neighborhood *pNbrhd = NULL);

      pRPL::Neighborhood* getNbrhd();
      const pRPL::Neighborhood* getNbrhd() const;

      void clearDataSettings();

      /* Option information */
      bool onlyUpdtCtrCell() const;
      bool needExchange() const;
      bool edgesFirst() const;

      /* Collective Evaluation  */
      pRPL::EvaluateReturn evalBR(const pRPL::CoordBR &br);
      pRPL::EvaluateReturn evalRandomly(const pRPL::CoordBR &br);
      pRPL::EvaluateReturn evalSelected(const pRPL::CoordBR &br,
                                        const pRPL::LongVect &vlclIdxs);

      /* ------------------------------------------------- *
       *             User-defined Processing               *
       * ------------------------------------------------- */

      /* After setting the Cellspaces for the Transition,
       * if the Cellspaces are SubCellspaces, the global ID of the SubCellspaces
       * will be given to this function.
       * The Transition may use this information to set other user-defined options
       * before running.
       * */
      virtual bool afterSetCellspaces(int subCellspcGlbIdx = pRPL::ERROR_ID);

      /* After setting the Neighborhood for the Transition,
       * the Transition may set other user-defined options before running
       * */
      virtual bool afterSetNbrhd();

      /* Final check before running */
      virtual bool check() const;

      /* Evaluate a specific Cell */
      virtual pRPL::EvaluateReturn evaluate(const pRPL::CellCoord &coord);

      /* Calculate the workload */
      virtual double workload(const pRPL::CoordBR &workBR) const;

    protected:
      string _primeLyrName;
      vector<string> _vInLyrNames;
      vector<string> _vOutLyrNames;
      map<string, pRPL::Cellspace*> _mpCellspcs;

      string _nbrhdName;
      pRPL::Neighborhood *_pNbrhd;

      bool _onlyUpdtCtrCell;
      bool _needExchange;
      bool _edgesFirst;
  };
};

#endif
