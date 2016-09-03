#include "prpl-transition.h"
#include <cmath>
#include <functional>

using namespace std;

namespace {
  static const double PI = 3.14159265;
};

class AspectTransition: public pRPL::Transition {
  public:
    AspectTransition();
    ~AspectTransition();

    void scale(float h2vScale);
    const float& scale() const;

    bool check() const;
    bool afterSetCellspaces(int subCellspcGlbIdx = pRPL::ERROR_ID);
    pRPL::EvaluateReturn evaluate(const pRPL::CellCoord &coord);

  protected:
    pRPL::Cellspace *_pDemCellspc;
    pRPL::Cellspace *_pSlpCellspc;
    pRPL::Cellspace *_pAspCellspc;
    int _demNoData;
    float _slpNoData;
    float _aspNoData;

    float _scale;
    float _cellWidth;
    float _cellHight;
};
