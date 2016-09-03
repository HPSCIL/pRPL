#ifndef PGTIOL_BASICTYPES_H        
#define PGTIOL_BASICTYPES_H        

#include <string>
#include <gdal_priv.h>
#include <gdal.h>

using namespace std;

#define MAXLEN 200                        
#define DATAOFFSET 8
#define TILEWIDTH 128           //default  for tilewidth

//geotiff type
enum GTiffType {
  TILETYPE,
  STRIPTYPE
};

struct Rectangle{
    Rectangle()
    :pXstart(0),pYstart(0),pnumCols(0),pnumRows(0){};
    Rectangle(int64_t Xstart,
        int64_t Ystart,
        int64_t numCols,
        int64_t numRows)
    :pXstart(Xstart),pYstart(Ystart),pnumCols(numCols),pnumRows(numRows) {};
    int64_t pXstart;
    int64_t pYstart;
    int64_t pnumCols;
    int64_t pnumRows;
};

struct Coordinate {
    Coordinate()
      :x(0.0), y(0.0) {}
    Coordinate(double xx, double yy)
      :x(xx), y(yy) {}
    bool operator==(const Coordinate &c) {
      return (x == c.x && y == c.y);
    }
    bool operator!=(const Coordinate &c) { return !(*this == c); }
    double x;
    double y;
};

struct Area {
    Area()
      :ul(), lr() {}
    Area(double ulx,
         double uly,
         double lrx,
         double lry)
      :ul(ulx, uly), lr(lrx, lry){}
    Coordinate ul;
    Coordinate lr;
};

//tiff metadata
struct geoTiffMetaData{
    double geoTransform[6];
    string projectionRef;
    int64_t x_size;
    int64_t y_size;
    int band_count;
    double no_data;
    GDALDataType band_type;
    int band_type_size;
    int64_t first_strip_offset;  //
    int64_t *tile_offsets;
    long blocks_x_size;    //tile or strip 
    long blocks_y_size;   //tile or strip 
    int64_t tiles_across;
    int64_t tiles_down;
    GDALColorTable ct;//colortable
    unsigned short sampleFormat;    
    unsigned short dataSizeObj;     
    unsigned short dataSizeFileIn; 
    unsigned short isCompressed;
    GTiffType tileOrStrip;
   
};

#endif 
