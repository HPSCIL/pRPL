#ifndef  PGTIOL_GTIFFDATA_H
#define  PGTIOL_GTIFFDATA_H

#include <gdal.h>
#include <gdal_priv.h>
#include <mpi.h>
#include <vector>
#include "math.h"
#include "pgtiol-basicTypes.h"
#include "pgtiol-gTiffMetaData.h"

using namespace std;

class gTiffData {
  public:
    gTiffData(void);
    ~gTiffData(void);
    bool setInputFileName(char*fileName);
    bool setOutputFileName(char*fileName);
    bool readDataFromGDAL(long xstart,long ystart, long numRows,long numCols,void* dest,GDALDataType type,int iband);
    bool readData(long xstart,
                  long ystart,
                  long numRows,
                  long numCols,
                  void* dest,
                  GDALDataType type,
                  int iband,
                  gTiffMetaData*iMeataData);

    void initWrite();
    bool writeData(long xstart,
                   long ystart,
                   long numRows,
                   long numCols,
                   void* source,
                   GDALDataType type,
                   gTiffMetaData*iMeataData,
				   long nLineSpace);

  private:
    char inputgTiff[MAXLEN];
    char outputgTiff[MAXLEN];
    int rank;

    bool write_area(gTiffMetaData*iMeataData,
                    void *data,
                    double ul_x,
                    double ul_y,
                    double lr_x,
                    double lr_y,
					long nLineSpace);

    Area calculate_tile_intersection(gTiffMetaData*iMeataData, Area subset);
    bool fill_stack(std::vector<Area> *write_stack,
                    Area old_area,
                    Area written_subset);

    bool write_subset(MPI_File&fh,gTiffMetaData*iMeataData,
                      void *data,
                      int64_t buffer_ul_x,
                      int64_t buffer_ul_y,
                      int64_t buffer_lr_x,
                      int64_t buffer_lr_y,
                      int64_t write_ul_x,
                      int64_t write_ul_y,
                      int64_t write_lr_x,
                      int64_t write_lr_y,
					  long nLineSpace);

    int64_t calculate_file_offset(gTiffMetaData*iMeataData,
                                  const int64_t raster_x,
                                  const int64_t raster_y);

    bool readStripData(long xstart,
                       long ystart,
                       long numRows,
                       long numCols,
                       void* dest,
                       gTiffMetaData*iMeataData);

    void getStripRowData(void* dest,
                         MPI_File&fh,
                         long currentStripFirstRowIndex,
                         long currentStripLastRowIndex,
                         long currentStripIndex,
                         long readCols,
                         uint32_t stripWidth,
                         uint32_t stripLength,
                         unsigned short sampleFormat,
                         unsigned short dataSizeFileIn,
                         double NoData,
                         GDALDataType dataType,
                         Rectangle rectangle);

    template <class Soutype,class Destype>
    void readStripByRow(MPI_File&fh,
                        Soutype type1,
                        Destype type2,
                        long currentStripIndex,
                        long currentRowIndex,
                        long readCols,
                        uint32_t stripWidth,
                        uint32_t stripLength,
                        unsigned short dataSizeFileIn,
                        Rectangle rectangle,
                        double NoData,
                        Soutype*tempbuff,
                        void*dest);
						
						
	bool  readTileData(long xstart,
                       long ystart,
                       long numRows,
                       long numCols,
                       void* dest,
                       gTiffMetaData*iMeataData);

    template <class memType>						  
    memType*getMem(memType type,long memSize);

    template <class Soutype,class Destype>
    void readTileDataByRow(MPI_File&fh,
                           Soutype type1,
                           Destype type2,
                           long tileCols,
                           uint32_t destOffset,
                           long tileWidth,
                           unsigned short dataSizeFileIn,
                           double NoData,
                           Soutype*tempbuff,
                           long readBytes,
                           void*dest);

    bool readTileDataByRows(long xstart,
                            long ystart,
                            long numRows,
                            long numCols,
                            void* dest,
                            gTiffMetaData*iMeataData);
							
	bool readTileDataByCols(long xstart, 
						    long ystart,
						    long numRows,
						    long numCols,
						    void* dest,
						    gTiffMetaData*iMeataData);
    template <class Soutype,class Destype>
    bool readTileDataByCol( Soutype type1,
							Destype type2,
							uint32_t ReadstartX,
							uint32_t ReadendX,
							uint32_t tileCols,
							uint32_t destOffset,
							unsigned short dataSizeFileIn,
							Soutype *tempDataRow,
							MPI_File&fh,
							void*dest);

};

#endif

