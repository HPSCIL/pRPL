#include "pgtiol-gTiffData.h"

gTiffData::gTiffData(void) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}

gTiffData::~gTiffData(void) {}

bool gTiffData::setInputFileName(char*fileName) {
  if(fileName == NULL){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: input filename is null (" << fileName << ")" << endl;
    return false;
  }

  strcpy(inputgTiff, fileName);
  return true;
}

bool gTiffData::setOutputFileName(char*fileName) {
  if(fileName==NULL){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: output filename is null (" << fileName << ")" << endl;
    return false;
  }

  strcpy(outputgTiff, fileName);
  return true;
}

bool gTiffData::readDataFromGDAL(long xstart, 
		                             long ystart,
		                             long numRows,
		                             long numCols,
		                             void* dest,
		                             GDALDataType type,
		                             int iband) {
  GDALAllRegister();
  GDALDataset *input_raster =static_cast<GDALDataset*>(GDALOpen(inputgTiff, GA_ReadOnly));
  if (input_raster == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << "Error opening input raster"<<endl;
    return false;
  }
  GDALRasterBand *pBand = input_raster->GetRasterBand(iband);
  if (pBand->RasterIO(GF_Read, xstart, ystart, numCols, numRows,
                      dest, numCols, numRows, type, 0, 0) != CE_None){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << "error  happen on read geotiff data"<<endl;
    return false;
  }
  GDALClose(input_raster);
  return true;
}

bool gTiffData::readData(long xstart, 
                         long ystart,
                         long numRows,
                         long numCols,
                         void* dest,
                         GDALDataType type,
                         int iband,
                         gTiffMetaData*iMeataData) {
  // cout<<iMeataData->tiffMetaData->blocks_x_size<<"read"<<iMeataData->tiffMetaData->blocks_y_size<<endl;
  if(iMeataData->getMetaData()->isCompressed == 1) {
    if(iMeataData->getMetaData()->tileOrStrip == 1) { 
      if(!readStripData(xstart, ystart, numRows, numCols,dest,iMeataData)){
        cerr << __FILE__ << " " << __FUNCTION__ \
            << " Error occurred when reading geotiff data" << endl;
        return false;
      }
    }
    else { 
      if(!readTileData(xstart, ystart, numRows, numCols,dest,iMeataData)){
        cerr << __FILE__ << " " << __FUNCTION__ \
            << " Error occurred when reading geotiff data"<<endl;
        return false;
      }
    }
	/*if(!readDataFromGDAL(xstart,ystart, numRows,numCols,dest,type,iband)){
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error occurred when reading geotiff data"<<endl;
      return false;
	}*/
  }
  else { //compressed
    if(!readDataFromGDAL(xstart,ystart, numRows,numCols,dest,type,iband)){
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error occurred when reading geotiff data"<<endl;
      return false;
    }
  }
  return true;
}

void gTiffData::initWrite()
{
  return;
}

/*bool gTiffData::writeData(long xstart, 
				          long ystart,
				          long numRows,
				          long numCols,
				          void* source,
						  GDALDataType type,
						  gTiffMetaData*iMeataData)
{
	cout<<"sssssssssaaaaaaaaaaaaaaaaa"<<outputgTiff<<endl;
	MPI_File fh;
    int	file_error = MPI_File_open( MPI_COMM_SELF, outputgTiff, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	if(file_error != MPI_SUCCESS){
		MPI_Abort(MPI_COMM_WORLD,22);
		return false;
	}
	int _rank=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&_rank);

	if(_rank==1){
	     for(int i=0;i<10;i++){
	        cout<<"rankwrite"<<_rank<<((float*)source)[i];
	    }
	    cout<<endl;
     }


    cout<<_rank<<"iMeataData->tiffMetaData->x_size"<<iMeataData->tiffMetaData->x_size    
	           <<"iMeataData->tiffMetaData->y_size"<<iMeataData->tiffMetaData->y_size    
	           <<"iMeataData->tiffMetaData->band_type"<<iMeataData->tiffMetaData->band_type
	           <<"iMeataData->tiffMetaData->x_size"<<iMeataData->tiffMetaData->dataSizeObj
			   <<"xstart" <<xstart
			   <<"ystart" <<ystart
			   <<"numRows"<<numRows
			   <<"numCols"<<numCols
			   <<"DATAOFFSET"<<DATAOFFSET<<endl;


	long partOffset = 0;
	long domainX =numCols;
	MPI_Status status;
	MPI_Offset mpiOffset;
	long totalX =iMeataData->tiffMetaData->x_size; 
	long totalY =iMeataData->tiffMetaData->y_size;
	long XstartNum=(totalX*ystart) + xstart;
	long next = partOffset;
	GDALDataType dataType  = iMeataData->tiffMetaData->band_type;
	unsigned short dataSizeObj = iMeataData->tiffMetaData->dataSizeObj;
	if(dataType == GDT_Int16) {        
		mpiOffset = DATAOFFSET + (dataSizeObj*(XstartNum));
		MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
		for(long i=0; i<numRows; i++)
		{       
			MPI_File_write(fh, (void*)(((short*)source)+next), numCols, MPI_SHORT, &status);
			mpiOffset += (totalX*dataSizeObj);
			MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
			next+=domainX;
		}
	}
	else if(dataType == GDT_Int32 ) 
	{
		mpiOffset = DATAOFFSET + (dataSizeObj*(XstartNum));
		MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
		int32_t* tempDataRow = new int32_t[numCols];
		for(long i=0; i<numRows; i++) 
		{ 
			for(long k = 0; k < numCols; k++)
			{
				tempDataRow[k]=(int32_t) (*((long*)source+next+k));
		    }
			MPI_File_write( fh, (void*)tempDataRow, numCols*dataSizeObj, MPI_BYTE, &status);
			mpiOffset += (totalX*dataSizeObj);
			MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
			next+=domainX;
		}
		delete tempDataRow;
	}
	else if( dataType == GDT_Float32 )
	{
		cout<<"float"<<endl;
		mpiOffset = DATAOFFSET + (dataSizeObj*(XstartNum));
		MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
		for(long i=0; i<numRows; i++) 
		{   

			MPI_File_write(fh, (void*)(((float*)source)+next), numCols, MPI_FLOAT, &status);
			mpiOffset += (totalX*dataSizeObj);
			MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
			next+=domainX;
		}
	}
	else if(dataType==GDT_UInt16)
	{
		mpiOffset = DATAOFFSET + (dataSizeObj*(XstartNum));
		MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
		for(long i=0; i<numRows; i++) 
		{
			MPI_File_write( fh, (void*)(((unsigned short*)source)+next), numCols, MPI_UNSIGNED_SHORT, &status);
			mpiOffset += (totalX*dataSizeObj);
			MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
			next+=domainX;
		}
	}
	else if(dataType==GDT_Int32)
	{
		mpiOffset = DATAOFFSET + (dataSizeObj*(XstartNum));
		MPI_File_seek(fh, mpiOffset, MPI_SEEK_SET);
		for(long i=0; i<numRows; i++) 
		{

			MPI_File_write(fh, (void*)(((int*)source)+next), numCols, MPI_INT, &status);
			mpiOffset += (totalX*dataSizeObj);
			MPI_File_seek( fh, mpiOffset, MPI_SEEK_SET);
			next+=domainX;
		}
	}
	MPI_File_close(&fh);
	return true;
}*/

bool gTiffData::writeData(long xstart, 
                          long ystart,
                          long numRows,
                          long numCols,
                          void* source,
                          GDALDataType type,
                          gTiffMetaData*iMeataData,
						  long nLineSpace) {
  Area write_area1(xstart,ystart,xstart+numCols-1,ystart+numRows-1);
  bool err = write_area(iMeataData,source,write_area1.ul.x,write_area1.ul.y,write_area1.lr.x,write_area1.lr.y,nLineSpace);
  if(err == false) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error occurred when writing geotiff data"<<endl;
    MPI_Abort(MPI_COMM_WORLD,22);
    return false;
  }
  return true;
}

bool gTiffData::write_area(gTiffMetaData*iMeataData,
                           void *data,
                           double ul_x,
                           double ul_y,
                           double lr_x,
                           double lr_y,
						   long nLineSpace) {
  
  std::vector<Area> write_stack;
  Area write_area;
  write_area.ul = Coordinate(ul_x, ul_y);
  write_area.lr = Coordinate(lr_x, lr_y);
  MPI_File fh;
  int rc = MPI_File_open(MPI_COMM_SELF,
                         outputgTiff,
                         MPI_MODE_RDWR,
                         MPI_INFO_NULL,
                         &fh);
  if (rc != MPI_SUCCESS) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: MPI_file unable to open file (" << outputgTiff << ")" << endl;
    return false;
  }

  write_stack.push_back(write_area);

  while (!write_stack.empty()) {
    Area top = write_stack.back();
    write_stack.pop_back();
    int rank = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Area subset = calculate_tile_intersection(iMeataData, top);
    fill_stack(&write_stack, top, subset);
    write_subset(fh,
	             iMeataData,
                 data,
                 ul_x,
                 ul_y,
                 lr_x,
                 lr_y,
                 subset.ul.x,
                 subset.ul.y,
                 subset.lr.x,
                 subset.lr.y,
				 nLineSpace);
  }
  MPI_File_close(&fh);
  return true;
}

Area gTiffData::calculate_tile_intersection(gTiffMetaData*iMeataData,
                                            Area subset) {
  const double tile_x_beginning = (static_cast<int64_t>(subset.ul.x)
                                   / iMeataData->getMetaData()->blocks_x_size)
                                   * iMeataData->getMetaData()->blocks_x_size;
  const double tile_y_beginning = (static_cast<int64_t>(subset.ul.y)
                                   / iMeataData->getMetaData()->blocks_y_size)
                                   * iMeataData->getMetaData()->blocks_y_size;
  const double tile_x_end = tile_x_beginning + iMeataData->getMetaData()->blocks_x_size - 1;
  const double tile_y_end = tile_y_beginning + iMeataData->getMetaData()->blocks_y_size - 1;
  const double subset_lr_x = min(tile_x_end, subset.lr.x);
  const double subset_lr_y = min(tile_y_end, subset.lr.y);

  return Area(subset.ul.x,
              subset.ul.y,
              subset_lr_x,
              subset_lr_y);
}

bool gTiffData::fill_stack(std::vector<Area> *write_stack,
                           Area old_area,
                           Area written_subset) {
  const double size_below = old_area.lr.y - written_subset.lr.y;
  const double size_right = old_area.lr.x - written_subset.lr.x;

  if (size_right > 0.0) {
    write_stack->push_back(Area(written_subset.lr.x + 1,
                           old_area.ul.y,
                           old_area.lr.x,
                           old_area.lr.y));
  }

  if (size_below > 0.0) {
    write_stack->push_back(Area(old_area.ul.x,
                           written_subset.lr.y + 1,
                           written_subset.lr.x,
                           old_area.lr.y));
  }
  return true;
}

bool gTiffData::write_subset(MPI_File&fh,
                             gTiffMetaData*iMeataData,
                             void *data,
                             int64_t buffer_ul_x,
                             int64_t buffer_ul_y,
                             int64_t buffer_lr_x,
                             int64_t buffer_lr_y,
                             int64_t write_ul_x,
                             int64_t write_ul_y,
                             int64_t write_lr_x,
                             int64_t write_lr_y,
							 long nLineSpace) {
  const double tile_x_beginning = (write_ul_x / iMeataData->getMetaData()->blocks_x_size)
                                  * iMeataData->getMetaData()->blocks_x_size;
  const double tile_x_end = tile_x_beginning + iMeataData->getMetaData()->blocks_x_size - 1;
  int count = 0;
  MPI_Status status;
  if (write_ul_x == tile_x_beginning
      && write_lr_x == tile_x_end) {
    count = ((write_lr_x - write_ul_x + 1) * (write_lr_y - write_ul_y + 1)
            * iMeataData->getMetaData()->band_type_size * iMeataData->getMetaData()->band_count);

    const int sub_row_size = (write_lr_x - write_ul_x + 1)
                            * iMeataData->getMetaData()->band_type_size
                            * iMeataData->getMetaData()->band_count;

    char *buffer = new(std::nothrow) char[count];
    for (int y = write_ul_y; y <= write_lr_y; ++y) {
      int pixel_offset = ((y - buffer_ul_y) * nLineSpace+ (write_ul_x - buffer_ul_x))
              * iMeataData->getMetaData()->band_type_size * iMeataData->getMetaData()->band_count;

      memcpy(buffer+((y-write_ul_y)*sub_row_size),
             static_cast<char*>(data)+pixel_offset,
             sub_row_size);
    }

    MPI_File_write_at(fh,
                      calculate_file_offset(iMeataData, write_ul_x, write_ul_y),
                      buffer,
                      count,
                      MPI_BYTE,
                      &status);
    delete[] buffer;
  }

  else {
    count = (write_lr_x - write_ul_x + 1)
            * iMeataData->getMetaData()->band_type_size
            * iMeataData->getMetaData()->band_count;

    for (int y = write_ul_y; y <= write_lr_y; ++y) {
      int pixel_offset = ((y - buffer_ul_y) * nLineSpace
          + (write_ul_x - buffer_ul_x))
              * iMeataData->getMetaData()->band_type_size * iMeataData->getMetaData()->band_count;
      MPI_File_write_at(fh,
                        calculate_file_offset(iMeataData, write_ul_x, y),
                        (static_cast<char*>(data)) + pixel_offset,
                        count,
                        MPI_BYTE,
                        &status);
    }
  }
  return true;
}

int64_t gTiffData::calculate_file_offset(gTiffMetaData*iMeataData,
                                         const int64_t raster_x,
                                         const int64_t raster_y) {
  const int64_t tile_x = raster_x % iMeataData->getMetaData()->blocks_x_size;
  const int64_t tile_y = raster_y % iMeataData->getMetaData()->blocks_y_size;
  const int64_t offset_into_tile = (tile_x + (tile_y * iMeataData->getMetaData()->blocks_x_size))
          * iMeataData->getMetaData()->band_type_size
          * iMeataData->getMetaData()->band_count;
  const int64_t tile_index = (raster_x / iMeataData->getMetaData()->blocks_x_size)
          + (raster_y / iMeataData->getMetaData()->blocks_y_size) * iMeataData->getMetaData()->tiles_across;

  return iMeataData->getMetaData()->tile_offsets[tile_index] + offset_into_tile;
}

bool gTiffData::readStripData(long xstart,
                              long ystart,
                              long numRows,
                              long numCols,
                              void* dest,
                              gTiffMetaData*iMeataData) {
  MPI_Offset mpiOffset;
  Rectangle rectangle(xstart,ystart, numCols,numRows);

  long totalX = iMeataData->getRasterXSize();
  long totalY = iMeataData->getRasterYSize();
  long tileWidth=iMeataData->getTileWidth();
  long tileLength=iMeataData->getTileLength();
  unsigned short sampleFormat = iMeataData->getSampleFormat();
  unsigned short dataSizeFileIn =iMeataData->getDataSizeFileIn();
  GDALDataType dataType = iMeataData->getDataType();
  double NoData = iMeataData->getNoDataValue();

  long firstStripIndex = rectangle.pYstart/tileLength;;
  long lastStripIndex =  (rectangle.pYstart + rectangle.pnumRows-1)/tileLength;
  long firstStripFirstRowIndex =   rectangle.pYstart - ( firstStripIndex * tileLength);
  long lastStripLastRowIndex = rectangle.pYstart + rectangle.pnumRows-1 - ( lastStripIndex * tileLength);
  long currentStripFirstRowIndex=0;
  long currentStripLastRowIndex=0;

  /*	cout<<"tileWidth"<<tileWidth
	    <<"tileLength"<<tileLength
		<<"sampleFormat"<<sampleFormat
		<<"dataSizeFileIn"<<dataSizeFileIn
		<<"NoData"<<NoData
		<<"firstStripIndex"<<firstStripIndex
		<<"lastStripIndex"<<lastStripIndex
		<<"firstStripFirstRowIndex"<<firstStripFirstRowIndex
		<<"lastStripLastRowIndex"<<lastStripLastRowIndex
		<<"rectangle.pXstart"<<rectangle.pXstart
		<<"rectangle.pYstart"<<rectangle.pYstart
		<<"rectangle.pnumRows"<<rectangle.pnumRows
		<<"rectangle.pnumCols"<<rectangle.pnumCols<<endl;*/

  MPI_File fh;
  int rc = MPI_File_open(MPI_COMM_SELF,
                         inputgTiff,
                         MPI_MODE_RDWR,
                         MPI_INFO_NULL,
                         &fh);

  if (rc != MPI_SUCCESS) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: MPI_file unable to open file (" << inputgTiff << ")" << endl;
    return false;
  }

  for(long i = firstStripIndex; i <= lastStripIndex; ++i ){
    if( i > firstStripIndex){
      currentStripFirstRowIndex=0;
    }
    else{
      currentStripFirstRowIndex=firstStripFirstRowIndex;
    }
    if(i < lastStripIndex){
      currentStripLastRowIndex=tileLength-1;
    }
    else{
      currentStripLastRowIndex = lastStripLastRowIndex;
    }
    if(currentStripFirstRowIndex <= currentStripLastRowIndex){
      mpiOffset = iMeataData->getMetaData()->tile_offsets[i];
      MPI_File_seek( fh, mpiOffset, MPI_SEEK_SET);
      mpiOffset = dataSizeFileIn * tileWidth * currentStripFirstRowIndex;
      MPI_File_seek( fh, mpiOffset, MPI_SEEK_CUR);
      getStripRowData(dest,fh,currentStripFirstRowIndex,currentStripLastRowIndex,i,
          numCols,tileWidth,tileLength,sampleFormat,dataSizeFileIn,NoData,dataType,rectangle);
    }
  }
  MPI_File_close(&fh);
  return true;
}

void gTiffData::getStripRowData(void* dest,
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
                                Rectangle rectangle) {
  for(long j = currentStripFirstRowIndex; j <= currentStripLastRowIndex; ++j){
    if((sampleFormat == 1) && (dataSizeFileIn == 1)){
      uint8_t type =0;
      uint8_t *tempDataRow= getMem(type,stripWidth);

      if(dataType == GDT_Byte){
        uint8_t type1=0;
        unsigned char type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Int16){
        uint8_t type1=0;
        short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_UInt16){
        uint8_t type1=0;
        unsigned short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_Int32){
        uint8_t type1=0;
        int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_UInt32){
        uint8_t type1=0;
        unsigned int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Float32){
        uint8_t type1=0;
        float type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	  else if(dataType == GDT_Float64){
		uint8_t type1=0;
        double type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest); 
	  }
      delete []tempDataRow;
    }
    else if ((sampleFormat == 1) && (dataSizeFileIn == 2)){
      uint16_t type =0;
      uint16_t *tempDataRow= getMem(type,stripWidth);
      if(dataType == GDT_Byte){
        uint16_t type1=0;
        unsigned char type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Int16){
        uint16_t type1=0;
        short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_UInt16){
        uint16_t type1=0;
        unsigned short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_Int32){
        uint16_t type1=0;
        int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_UInt32){
        uint16_t type1=0;
        unsigned int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Float32){
        uint16_t type1=0;
        float type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	  else if(dataType == GDT_Float64){
        uint16_t type1=0;
        double type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      delete [] tempDataRow;
    }
    else if ((sampleFormat == 1) && (dataSizeFileIn == 4)) {
      uint32_t type =0;
      uint32_t *tempDataRow= getMem(type,stripWidth);
      if(dataType == GDT_Byte){
        uint32_t type1=0;
        unsigned char type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Int16){
        uint32_t type1=0;
        short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_UInt16){
        uint32_t type1=0;
        unsigned short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_Int32){
        uint32_t type1=0;
        int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_UInt32){
        uint32_t type1=0;
        unsigned int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Float32){
        uint32_t type1=0;
        float type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	   else if(dataType == GDT_Float64){
        uint32_t type1=0;
        double type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      delete []tempDataRow;
    }
    else if ((sampleFormat == 2) && (dataSizeFileIn == 1)){
      int8_t type =0;
      int8_t *tempDataRow= getMem(type,stripWidth);

      if(dataType == GDT_Byte){
        int8_t type1=0;
        unsigned char type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Int16){
        int8_t type1=0;
        short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_UInt16){
        int8_t type1=0;
        unsigned short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_Int32){
        int8_t type1=0;
        int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_UInt32){
        int8_t type1=0;
        unsigned int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Float32){
        int8_t type1=0;
        float type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	  else if(dataType == GDT_Float64){
        int8_t type1=0;
        double type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	  
      delete []tempDataRow;
    }
    else if ((sampleFormat == 2) && (dataSizeFileIn == 2)) {
      int16_t type =0;
      int16_t *tempDataRow= getMem(type,stripWidth);
      if(dataType == GDT_Byte){
        int16_t type1=0;
        unsigned char type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Int16){
        int16_t type1=0;
        short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_UInt16){
        int16_t type1=0;
        unsigned short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_Int32){
        int16_t type1=0;
        int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_UInt32){
        int16_t type1=0;
        unsigned int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Float32){
        int16_t type1=0;
        float type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	  else if(dataType == GDT_Float64){
        int16_t type1=0;
        double type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
    }
    else if ((sampleFormat == 2) && (dataSizeFileIn == 4)){
      int32_t type =0;
      int32_t *tempDataRow= getMem(type,stripWidth);

      if(dataType == GDT_Byte){
        int32_t type1=0;
        unsigned char type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Int16){
        int32_t type1=0;
        short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_UInt16){
        int32_t type1=0;
        unsigned short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_Int32){
        int32_t type1=0;
        int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_UInt32){
        int32_t type1=0;
        unsigned int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Float32){
        int32_t type1=0;
        float type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	  else if(dataType == GDT_Float64){
        int32_t type1=0;
        double type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	  
      delete []tempDataRow;
    }
    else if ((sampleFormat == 3) && (dataSizeFileIn == 4)) {
      float type =0;
      float *tempDataRow= getMem(type,stripWidth);

      if(dataType == GDT_Byte){
        float type1=0;
        unsigned char type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Int16){
        float type1=0;
        short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_UInt16){
        float type1=0;
        unsigned short type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_Int32){
        float type1=0;
        int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType ==GDT_UInt32){
        float type1=0;
        unsigned int type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      else if(dataType == GDT_Float32){
        float type1=0;
        float type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
	  else if(dataType == GDT_Float64){
        float type1=0;
        double type2=0;
        readStripByRow(fh,type1,type2,currentStripIndex,j,readCols,stripWidth,stripLength,dataSizeFileIn,rectangle,NoData,tempDataRow,dest);
      }
      delete []tempDataRow;
    }
    else{
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error: Unsupported TIFF file type. sampleFormat = "
          << sampleFormat << "; dataSizeFileIn = " << dataSizeFileIn << endl;
      MPI_Abort(MPI_COMM_WORLD,-1);
    }
  }
}

template <class Soutype,class Destype>
void gTiffData::readStripByRow(MPI_File&fh,
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
                               void*dest) {
    MPI_Status status;
    MPI_File_read(fh, tempbuff, dataSizeFileIn * stripWidth, MPI_BYTE, &status);
    for(long k =rectangle.pXstart; k<rectangle.pXstart+readCols; k++){
      ((Destype*)dest)[(currentRowIndex+(currentStripIndex*stripLength)-rectangle.pYstart)*readCols+(k-rectangle.pXstart)]=(Destype)(tempbuff[k]);
    }
}

bool  gTiffData::readTileData(long xstart,
                              long ystart,
                              long numRows,
                              long numCols,
                              void* dest,
                              gTiffMetaData*iMeataData){
						  
		if(numCols==iMeataData->getRasterXSize()){
			if(!readTileDataByRows(xstart, ystart, numRows, numCols,dest,iMeataData)){
                 cerr << __FILE__ << " " << __FUNCTION__ \
                      << " Error occurred when reading geotiff data"<<endl;
                 return false;
            }
		}else{
			if(!readTileDataByCols(xstart, ystart, numRows, numCols,dest,iMeataData)){
                 cerr << __FILE__ << " " << __FUNCTION__ \
                      << " Error occurred when reading geotiff data"<<endl;
                 return false;
            }
		}		
		return true;							  
}


template <class memType>
memType*gTiffData::getMem(memType type,long memSize){
    return new memType[memSize];
}

bool gTiffData::readTileDataByRows(long xstart, 
                                   long ystart,
                                   long numRows,
                                   long numCols,
                                   void* dest,
                                   gTiffMetaData*iMeataData) {
									   
									   									   
  long totalX = iMeataData->getRasterXSize();
  long totalY = iMeataData->getRasterYSize();
  long tileWidth=iMeataData->getTileWidth();
  long tileLength=iMeataData->getTileLength();
  unsigned short sampleFormat = iMeataData->getSampleFormat();
  unsigned short dataSizeFileIn =iMeataData->getDataSizeFileIn();
  GDALDataType dataType = iMeataData->getDataType();
  double NoData = iMeataData->getNoDataValue();
  
  MPI_Offset mpiOffset;
  uint32_t tileStart=0;
  uint32_t tilesAcross=(totalX-1)/tileWidth+1;
  uint32_t tileEnd=0;
  uint32_t rowInTile=0;
  uint32_t tileCols=0;
  uint32_t destOffset1=0;
  double timespendonmpiread=0;

  MPI_File fh;
  int rc = MPI_File_open(MPI_COMM_SELF,
                         inputgTiff,
                         MPI_MODE_RDONLY,
                         MPI_INFO_NULL,
                         &fh);

  if (rc != MPI_SUCCESS) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: MPI_file unable to open file (" << inputgTiff << ")" << endl;
    return false;
  }

  for(uint32_t y=ystart; y<ystart+numRows; y++) {
    tileStart=y/tileLength*tilesAcross;
    if(tileStart > 0){
      tileStart = tileStart;
    }
    tileEnd=tileStart+tilesAcross-1;
    rowInTile=y-(tileStart/tilesAcross)*tileLength;
    if(tileStart <= tileEnd) {
      //MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
      for(uint32_t tile=tileStart; tile<=tileEnd;tile++)
      {
        tileCols = tileWidth;
        if(tile == tileEnd)tileCols=totalX-(tilesAcross-1)*tileWidth;
        mpiOffset=iMeataData->getMetaData()->tile_offsets[tile]+rowInTile*tileWidth*dataSizeFileIn;

        long readBytes =dataSizeFileIn * tileCols;

        MPI_File_seek(fh,mpiOffset,MPI_SEEK_SET);
        destOffset1=(tile-tileStart)*tileWidth+(y-ystart)*numCols;
        if((sampleFormat == 1) && (dataSizeFileIn == 1)){
          uint8_t type =0;
          uint8_t *tempDataRow= getMem(type,tileWidth);
          if(dataType == GDT_Byte){
            uint8_t type1=0;
            unsigned char type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Int16){
            uint8_t type1=0;
            short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Float32){
            uint8_t type1=0;
            float type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_UInt16){
            uint8_t type1=0;
            unsigned short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_Int32){
            uint8_t type1=0;
            int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_UInt32){
            uint8_t type1=0;
            unsigned int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
		  else if(dataType ==GDT_Float64){
            uint8_t type1=0;
            double type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          delete []tempDataRow;
        }
        else if ((sampleFormat == 1) && (dataSizeFileIn == 2)) {
          uint16_t type =0;
          uint16_t *tempDataRow= getMem(type,tileWidth);
          if(dataType == GDT_Byte){
            uint16_t type1=0;
            unsigned char type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Int16){
            uint16_t type1=0;
            short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Float32){
            uint16_t type1=0;
            float type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);

          }
          else if(dataType == GDT_UInt16){
            uint16_t type1=0;
            unsigned short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_Int32){
            uint16_t type1=0;
            int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_UInt32){
            uint16_t type1=0;
            unsigned int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
		  else if(dataType ==GDT_Float64){
            uint16_t type1=0;
            double type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          delete []tempDataRow;
        }
        else if ((sampleFormat == 1) && (dataSizeFileIn == 4)){
          uint32_t type =0;
          uint32_t *tempDataRow= getMem(type,tileWidth);
          if(dataType == GDT_Byte){
            uint32_t type1=0;
            unsigned char type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Int16){
            uint32_t type1=0;
            short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Float32){
            uint32_t type1=0;
            float type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);

          }
          else if(dataType == GDT_UInt16){
            uint32_t type1=0;
            unsigned short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_Int32){
            uint32_t type1=0;
            int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_UInt32){
            uint32_t type1=0;
            unsigned int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
		  else if(dataType ==GDT_Float64){
            uint32_t type1=0;
            double type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          delete []tempDataRow;
        }
        else if ((sampleFormat == 2) && (dataSizeFileIn == 1)) {
          int8_t type =0;
          int8_t *tempDataRow= getMem(type,tileWidth);
          if(dataType == GDT_Byte){
            int8_t type1=0;
            unsigned char type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Int16){
            int8_t type1=0;
            short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Float32){
            int8_t type1=0;
            float type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);

          }
          else if(dataType == GDT_UInt16){
            int8_t type1=0;
            unsigned short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_Int32){
            int8_t type1=0;
            int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_UInt32){
            int8_t type1=0;
            unsigned int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
		  else if(dataType ==GDT_Float64){
            int8_t type1=0;
            double type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          delete []tempDataRow;
        }
        else if ((sampleFormat == 2) && (dataSizeFileIn == 2)) {
          int16_t type =0;
          int16_t *tempDataRow= getMem(type,tileWidth);
          if(dataType == GDT_Byte){
            int16_t type1=0;
            unsigned char type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Int16){
            int16_t type1=0;
            short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Float32){
            int16_t type1=0;
            float type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);

          }
          else if(dataType == GDT_UInt16){
            int16_t type1=0;
            unsigned short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_Int32){
            int16_t type1=0;
            int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_UInt32){
            int16_t type1=0;
            unsigned int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
		  else if(dataType ==GDT_Float64){
            int16_t type1=0;
            double type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          delete []tempDataRow;
        }
        else if ((sampleFormat == 2) && (dataSizeFileIn == 4)) {
          int32_t type =0;
          int32_t *tempDataRow= getMem(type,tileWidth);
          if(dataType == GDT_Byte){
            int32_t type1=0;
            unsigned char type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Int16){
            int32_t type1=0;
            short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Float32){
            int32_t type1=0;
            float type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_UInt16){
            int32_t type1=0;
            unsigned short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_Int32){
            int32_t type1=0;
            int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);

          }
          else if(dataType ==GDT_UInt32){
            int32_t type1=0;
            unsigned int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
		   else if(dataType ==GDT_Float64){
            int32_t type1=0;
            double type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          delete []tempDataRow;
        }
        else if(sampleFormat == 3 && dataSizeFileIn == 4){
          float type =0;
          float *tempDataRow= getMem(type,tileWidth);
          if(dataType == GDT_Byte){
            float type1=0;
            unsigned char type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Int16){
            float type1=0;
            short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType == GDT_Float32){
            float type1=0;
            float type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);

          }
          else if(dataType == GDT_UInt16){
            float type1=0;
            unsigned short type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_Int32){
            float type1=0;
            int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          else if(dataType ==GDT_UInt32){
            float type1=0;
            unsigned int type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
		  else if(dataType ==GDT_Float64){
            float type1=0;
            double type2=0;
            readTileDataByRow(fh,type1,type2,tileCols,destOffset1,tileWidth,dataSizeFileIn,NoData,tempDataRow,readBytes,dest);
          }
          delete []tempDataRow;
        }
        else{
          cerr << __FILE__ << " " << __FUNCTION__ \
              << " Error: Unsupported TIFF file type. sampleFormat = "
              << sampleFormat << "; dataSizeFileIn = " << dataSizeFileIn << endl;
          MPI_Abort(MPI_COMM_WORLD,-1);
        }
      }
    }
  }

  MPI_File_close(&fh);
  return true;
}

template <class Soutype,class Destype>
void gTiffData::readTileDataByRow(MPI_File&fh,
                                  Soutype type1,
                                  Destype type2,
                                  long tileCols,
                                  uint32_t destOffset,
                                  long tileWidth,
                                  unsigned short dataSizeFileIn,
                                  double NoData,
                                  Soutype*tempbuff,
                                  long readBytes,
                                  void*dest) {
  MPI_Status status;
  MPI_File_read(fh, tempbuff, readBytes, MPI_BYTE, &status);
  for(long k = 0; k < tileCols; k++){
    ((Destype*)dest)[destOffset+k]=(Destype)(tempbuff[k]);
  }
}



bool gTiffData::readTileDataByCols(long xstart, 
                                   long ystart,
                                   long numRows,
                                   long numCols,
                                   void* dest,
                                   gTiffMetaData*iMeataData)
{
	MPI_Offset mpiOffset;
	uint32_t tileStart=0;
	uint32_t tilesAcross=0;
	uint32_t tileEnd=0;
	uint32_t rowInTile=0;
	uint32_t tileCols=0;
	uint32_t destOffset=0;
  
	long totalX = iMeataData->getRasterXSize();
    long totalY = iMeataData->getRasterYSize();
    long tileWidth=iMeataData->getTileWidth();
    long tileLength=iMeataData->getTileLength();
    unsigned short sampleFormat = iMeataData->getSampleFormat();
    unsigned short dataSizeFileIn =iMeataData->getDataSizeFileIn();
    GDALDataType dataType = iMeataData->getDataType();
    double NoData = iMeataData->getNoDataValue();
	uint32_t domainX =numCols;  
	MPI_File fh;
    int rc = MPI_File_open(MPI_COMM_SELF,
                           inputgTiff,
                           MPI_MODE_RDONLY,
                           MPI_INFO_NULL,
                           &fh);

    if (rc != MPI_SUCCESS) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: MPI_file unable to open file (" << inputgTiff << ")" << endl;
    return false;
    }
	
	tilesAcross=(totalX-1)/tileWidth+1; 
	uint32_t  AcrossStart=xstart/tileWidth;
	uint32_t  AcrossEnd=(xstart+numCols-1)/tileWidth;
	uint32_t  OverNum1=xstart%tileWidth;
	uint32_t  OverNum2=(xstart+numCols)%tileWidth;
	
		
	for(uint32_t y=ystart; y<ystart+numRows; y++)
	{
		tileStart=y/tileLength*tilesAcross;  
		if(tileStart > 0){
			tileStart=tileStart;
		}

		tileEnd=tileStart+tilesAcross-1;
		rowInTile=y-(tileStart/tilesAcross)*tileLength;
		uint32_t   AcrossStartX=y/tileLength*tilesAcross+AcrossStart;
		uint32_t   AcrossEndX=y/tileLength*tilesAcross+AcrossEnd;
		uint32_t   ReadstartX=0;
		uint32_t   ReadendX=0;

		if(AcrossStartX <= AcrossEndX){       
			for(uint32_t tile=AcrossStartX; tile<=AcrossEndX;tile++)
			{        
				tileCols = tileWidth;
				if(tile == tileEnd)tileCols=totalX-(tilesAcross-1)*tileWidth;
				mpiOffset=iMeataData->getMetaData()->tile_offsets[tile]+rowInTile*tileWidth*dataSizeFileIn;
				MPI_File_seek(fh,mpiOffset,MPI_SEEK_SET);
				if(tile==AcrossStartX)
					ReadstartX=OverNum1;
				else
					ReadstartX=0;
				if(tile==AcrossEndX)
					ReadendX=OverNum2;
				else
					ReadendX=tileCols;
				if(tile==AcrossStartX){
					destOffset=(y-ystart)*domainX;
				}
				else
					destOffset=(tile-AcrossStartX-1)*tileWidth+(tileWidth-OverNum1)+(y-ystart)*domainX;
				
				if((sampleFormat == 1) && (dataSizeFileIn == 1)) {
					uint8_t type =0;
                    uint8_t *tempDataRow= getMem(type,tileWidth);
					
					if(dataType == GDT_Byte){
						uint8_t type1=0;
						unsigned char type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					
					else if(dataType == GDT_Int16){
						uint8_t type1=0;
						short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int32){
						uint8_t type1=0;
						int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Float32){
						uint8_t type1=0;
						float type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_UInt16){
						uint8_t type1=0;
						unsigned short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_UInt32){
						uint8_t type1=0;
						unsigned int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_Float64){
						uint8_t type1=0;
						double type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
                    delete []tempDataRow;
				}
				else if ((sampleFormat == 1) && (dataSizeFileIn == 2)) {
					uint16_t type =0;
                    uint16_t *tempDataRow= getMem(type,tileWidth);
					
					if(dataType == GDT_Byte){
						uint16_t type1=0;
						unsigned char type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int16){
						uint16_t type1=0;
						short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int32){
						uint16_t type1=0;
						int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Float32){
						uint16_t type1=0;
						float type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);

					}
					else if(dataType == GDT_UInt16){
						uint16_t type1=0;
						unsigned short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_UInt32){
						uint16_t type1=0;
						unsigned int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_Float64){
						uint16_t type1=0;
						double type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					delete []tempDataRow;
				}
				else if ((sampleFormat == 1) && (dataSizeFileIn == 4)) {
					uint32_t type =0;
                    uint32_t *tempDataRow= getMem(type,tileWidth);
					
					if(dataType == GDT_Byte){
						uint32_t type1=0;
						unsigned char type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int16){
						uint32_t type1=0;
						short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int32){
						uint32_t type1=0;
						int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Float32){
						uint32_t type1=0;
						float type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);

					}
					else if(dataType == GDT_UInt16){
						uint32_t type1=0;
						unsigned short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_UInt32){
						uint32_t type1=0;
						unsigned int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_Float64){
						uint32_t type1=0;
						double type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					delete []tempDataRow;
				}
				else if ((sampleFormat == 2) && (dataSizeFileIn == 1)) {
					int8_t type =0;
                    int8_t *tempDataRow= getMem(type,tileWidth);
					
					if(dataType == GDT_Byte){
						int8_t type1=0;
						unsigned char type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int16){
						int8_t type1=0;
						short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int32){
						int8_t type1=0;
						int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Float32){
						int8_t type1=0;
						float type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);

					}
					else if(dataType == GDT_UInt16){
						int8_t type1=0;
						unsigned short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_UInt32){
						int8_t type1=0;
						unsigned int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_Float64){
						int8_t type1=0;
						double type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					delete []tempDataRow;
				}
				else if ((sampleFormat == 2) && (dataSizeFileIn == 2)){
					int16_t type =0;
                    int16_t *tempDataRow= getMem(type,tileWidth);
					
					if(dataType == GDT_Byte){
						int16_t type1=0;
						unsigned char type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int16){
						int16_t type1=0;
						short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int32){
						int16_t type1=0;
						int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Float32){
						int16_t type1=0;
						float type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);

					}
					else if(dataType == GDT_UInt16){
						int16_t type1=0;
						unsigned short type2=0;
					    readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_UInt32){
						int16_t type1=0;
						unsigned int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_Float64){
						int16_t type1=0;
						double type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					delete []tempDataRow;
				}
				else if ((sampleFormat == 2) && (dataSizeFileIn == 4)) {
			        int32_t type =0;
                    int32_t *tempDataRow= getMem(type,tileWidth);
					
					if(dataType == GDT_Byte){
						int32_t type1=0;
						unsigned char type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int16){
						int32_t type1=0;
						short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int32){
						int32_t type1=0;
						int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Float32){
						int32_t type1=0;
						float type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);

					}
					else if(dataType == GDT_UInt16){
						int32_t type1=0;
						unsigned short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_UInt32){
						int32_t type1=0;
						unsigned int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_Float64){
						int32_t type1=0;
						double type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					delete []tempDataRow;
				}
				else if(sampleFormat == 3 && dataSizeFileIn == 4){
					float type =0;
                    float *tempDataRow= getMem(type,tileWidth);
					if(dataType == GDT_Byte){
						float type1=0;
						unsigned char type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int16){
						float type1=0;
						short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Int32){
						float type1=0;
						int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType == GDT_Float32){
						float type1=0;
						float type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);

					}
					else if(dataType == GDT_UInt16){
						float type1=0;
						unsigned short type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					else if(dataType ==GDT_UInt32){
						float type1=0;
						unsigned int type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}else if(dataType ==GDT_Float64){
						float type1=0;
						double type2=0;
						readTileDataByCol(type1,type2,ReadstartX,ReadendX,tileCols,destOffset,dataSizeFileIn,tempDataRow,fh,dest);
					}
					delete []tempDataRow;
				}
				else{
				     cerr << __FILE__ << " " << __FUNCTION__ \
			             << "Unsupported TIFF file type.  _sampleFormat = "
			             <<sampleFormat<<"_dataSizeFileIn = "<<dataSizeFileIn<<endl;
					 
					 return false;
				}
			}


		}
	}
	MPI_File_close(&fh);
    return true;
}
template <class Soutype,class Destype>
bool gTiffData::readTileDataByCol(Soutype type1,
									Destype type2,
									uint32_t ReadstartX,
									uint32_t ReadendX,
									uint32_t tileCols,
									uint32_t destOffset,
									unsigned short dataSizeFileIn,
									Soutype *tempDataRow,
									MPI_File&fh,
									void*dest)
{  
    
	MPI_Status status;
	MPI_File_read( fh, tempDataRow, dataSizeFileIn * tileCols, MPI_BYTE, &status);
	for(long k =ReadstartX; k <ReadendX; k++)
	{ 
        
		((Destype*)dest)[destOffset+k-ReadstartX]=(Destype)tempDataRow[k];
	}
	return true;
}






