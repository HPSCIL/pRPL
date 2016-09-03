#include "aspectTrans.h"

AspectTransition::
AspectTransition()
  :pRPL::Transition(),
   _pDemCellspc(NULL),
   _pSlpCellspc(NULL),
   _pAspCellspc(NULL),
   _demNoData(pRPL::DEFAULT_NODATA_INT),
   _slpNoData(pRPL::DEFAULT_NODATA_FLOAT),
   _aspNoData(pRPL::DEFAULT_NODATA_FLOAT),
   _scale(1.0),
   _cellWidth(1.0),
   _cellHight(1.0) {
  _needExchange = false;
  _edgesFirst = false;
}

AspectTransition::
~AspectTransition() {
  _pDemCellspc = NULL;
  _pSlpCellspc = NULL;
  _pAspCellspc = NULL;
}

void AspectTransition::
scale(float h2vScale) {
  _scale = h2vScale;
}

const float& AspectTransition::
scale() const {
  return _scale;
}

bool AspectTransition::
check() const {
  if(_mpCellspcs.size() != 3) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: THREE Cellspaces are needed" \
         << endl;
    return false;
  }

  const pRPL::Cellspace *pPrmCellspc = getCellspaceByLyrName(_primeLyrName);
  if(pPrmCellspc == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: cannot find the primary Cellspace" \
         << endl;
    return false;
  }
  const pRPL::SpaceDims &myDims = pPrmCellspc->info()->dims();

  map<string, pRPL::Cellspace *>::const_iterator itrCellspcMap = _mpCellspcs.begin();
  while(itrCellspcMap != _mpCellspcs.end()) {
    const pRPL::Cellspace *pCellspc = itrCellspcMap->second;
    if(pCellspc == NULL) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: NULL pointer to Cellspace (" \
           << itrCellspcMap->first << ")" \
           << endl;
      return false;
    }
    if(pCellspc->info()->dims() != myDims) {
      cerr << __FILE__ << " function:" << __FUNCTION__ \
           << " Error: Cellspace(" << itrCellspcMap->first \
           << ")'s dimensions (" << pCellspc->info()->dims() \
           << ") do NOT match with the primary Cellspace's dimensions (" \
           << myDims << ")" \
           << endl;
      return false;
    }
    itrCellspcMap++;
  } // end -- while(itrCellspcMap != _mpCellspcs.end())

  if(_pNbrhd == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: NULL pointer to the Neighborhood" \
         << endl;
    return false;
  }

  if(_scale <= 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid scale factor (" << _scale \
         << ")" << endl;
    return false;
  }

  if(_cellWidth <= 0 || _cellHight <= 0) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
         << " Error: invalid Cell size (" << _cellWidth << ", " << _cellHight \
         << ")" << endl;
    return false;
  }

  return true;
}

bool AspectTransition::
afterSetCellspaces(int subCellspcGlbIdx) {
  _pDemCellspc = getCellspaceByLyrName(_vInLyrNames[0]);
  if(_pDemCellspc == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: NULL pointer to input DEM Cellspace (" \
        << _vInLyrNames[0] << ")" \
        << endl;
    return false;
  }
  _demNoData = _pDemCellspc->info()->getNoDataValAs<int>();

  if(_pDemCellspc->info()->isGeoreferenced(false)) {
    const pRPL::GeoCoord& cellSize = _pDemCellspc->info()->georeference()->cellSize();
    _cellWidth = fabs(cellSize.x()) * _scale;
    _cellHight = fabs(cellSize.y()) * _scale;
  }

  _pSlpCellspc = getCellspaceByLyrName(_vOutLyrNames[0]);
  if(_pSlpCellspc == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: NULL pointer to output SLOPE Cellspace (" \
        << _vOutLyrNames[0] << ")" \
        << endl;
    return false;
  }
  _slpNoData = _pSlpCellspc->info()->getNoDataValAs<float>();

  _pAspCellspc = getCellspaceByLyrName(_vOutLyrNames[1]);
  if(_pAspCellspc == NULL) {
    cerr << __FILE__ << " function:" << __FUNCTION__ \
        << " Error: NULL pointer to output ASPECT Cellspace (" \
        << _vOutLyrNames[1] << ")" \
        << endl;
    return false;
  }
  _aspNoData = _pAspCellspc->info()->getNoDataValAs<float>();

  return true;
}

pRPL::EvaluateReturn AspectTransition::
evaluate(const pRPL::CellCoord &coord) {
  float slope = 0.0, aspect = 0.0;

  pRPL::IntVect vDemVals;
  if(!_pDemCellspc->nbrhdValsAs<int>(vDemVals, *_pNbrhd, coord, true, true)) {
    return pRPL::EVAL_FAILED;
  }

  if(std::find(vDemVals.begin(), vDemVals.end(), _demNoData) != vDemVals.end()) {
    slope = _slpNoData;
    aspect = _aspNoData;
  }
  else {
    float demDiff;
    demDiff = static_cast<float>((vDemVals[1] + vDemVals[4] + vDemVals[4] + vDemVals[6]) -
                                 (vDemVals[3] + vDemVals[5] + vDemVals[5] + vDemVals[8]));
    float az = demDiff / (8.0 * _cellWidth);
    demDiff = static_cast<float>((vDemVals[6] + vDemVals[7] + vDemVals[7] + vDemVals[8]) -
                                 (vDemVals[1] + vDemVals[2] + vDemVals[2] + vDemVals[3]));
    float bz = demDiff / (8.0 * _cellHight);
    slope = pow((pow(double(az), 2.0) + pow(double(bz), 2.0)), 0.5);
    slope = atan(slope) * 180.0 / PI;

    if(slope >= -pRPL::EPSINON &&
       slope <= pRPL::EPSINON) {
      aspect = -1.0;
    }
    else {
      /* Determine the quadrant & calculate aspect */
      if(az < - pRPL::EPSINON) {
        if(bz < -pRPL::EPSINON) {
          aspect = (PI / 2.0) - atan(fabs((bz / az)));
        }
        else {
          aspect = (PI / 2.0) + atan(fabs((bz / az)));
        }
      } /* End of if(az < 0.0) */
      else {
        if(az > pRPL::EPSINON) {
          if(bz < -pRPL::EPSINON) {
            aspect = (3.0 * PI / 2.0) + atan(fabs((bz / az)));
          }
          else {
            aspect = (3.0 * PI / 2.0) - atan(fabs((bz / az)));
          }
        } /* End of if(az > 0.0) */
        else {
          /* az equals zero */
          if (bz < -pRPL::EPSINON) {
            aspect = 0.0;
          }
          else {
            aspect = PI;
          }
        } /* End of if(az == 0.0) */
      } /* End of if(az >= 0.0) */

      /* Convert from radians to degrees */
      aspect = 360.0 * (aspect / (2.0 * PI));
      /* Correct to north */
      aspect += 180.0;
      if(aspect >= 360.0) {
        aspect -= 360.0;
      }
    } // end -- calculating aspect when slope > 0
  } // end -- calculating aspect
 
  if(!_pSlpCellspc->updateCellAs<float>(coord, slope, true) ||
     !_pAspCellspc->updateCellAs<float>(coord, aspect, true)) {
    return pRPL::EVAL_FAILED;
  }

  return pRPL::EVAL_SUCCEEDED;
}
