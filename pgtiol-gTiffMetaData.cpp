#include "pgtiol-gTiffMetaData.h"

gTiffMetaData::gTiffMetaData(void) {
  tiffMetaData = new geoTiffMetaData();
  tiffMetaData->x_size=0;
  tiffMetaData->y_size=0;
  tiffMetaData->band_count=0;
  tiffMetaData->no_data=0.0;
  tiffMetaData->band_type=GDT_Unknown;
  tiffMetaData->band_type_size=0;
  tiffMetaData->first_strip_offset=0;
  tiffMetaData->blocks_x_size=0;
  tiffMetaData->blocks_y_size=0;
  tiffMetaData->tiles_across=0;
  tiffMetaData->tiles_down=0;
  tiffMetaData->tileOrStrip = TILETYPE;
  tiffMetaData->isCompressed = 1;
}

gTiffMetaData::~gTiffMetaData(void) {
  if(tiffMetaData!=NULL){
    delete tiffMetaData;
    tiffMetaData = NULL;
  }
}

//--------------------read and write api----------------

bool gTiffMetaData::setInputFile(char*fileName) {
  if(fileName==NULL){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error：input filename is null (" <<fileName<<")"<<endl;
    return false;
  }
  strcpy(inputGTiffFile, fileName);
  return true;
}

bool gTiffMetaData::readMetaData() {
  // gdal
  GDALAllRegister();
  GDALDataset *rasterData =static_cast<GDALDataset*>(GDALOpen(inputGTiffFile, GA_ReadOnly));

  if(rasterData == NULL){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error：unable to open the input file (" << inputGTiffFile << ")" \
        <<endl;
    return false;
  }

  tiffMetaData->x_size = rasterData->GetRasterXSize();
  tiffMetaData->y_size = rasterData->GetRasterYSize();
  tiffMetaData->band_count = rasterData->GetRasterCount();
  tiffMetaData->no_data = rasterData->GetRasterBand(1)->GetNoDataValue();
  tiffMetaData->band_type = rasterData->GetRasterBand(1)->GetRasterDataType();
  tiffMetaData->band_type_size = GDALGetDataTypeSize(tiffMetaData->band_type)/8;
  tiffMetaData->first_strip_offset = -1;
  tiffMetaData->blocks_x_size = tiffMetaData->x_size;
  tiffMetaData->blocks_y_size = tiffMetaData->y_size;
  rasterData->GetGeoTransform(tiffMetaData->geoTransform);
  tiffMetaData->projectionRef = rasterData->GetProjectionRef();
  //GDALColorTable ct = *(rasterData->GetRasterBand(1)->GetColorTable());//
  GDALClose(rasterData);

  //libtiff
  TIFF *tiffds = TIFFOpen(inputGTiffFile, "r");
  if(tiffds ==NULL){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: unable to open the input file (" << inputGTiffFile << ")" \
        <<endl;
    return false;
  }

  int64_t *offset = NULL;
  int64_t tiles_per_image = -1;
  int64_t *tiff_offsets = NULL;

  //get dataSizeFileIn and sampleFormat
  unsigned short sampleFormat =0;
  unsigned short dataSizeFileIn=0;
  unsigned short isCompression=1;

  TIFFGetField(tiffds, TIFFTAG_SAMPLEFORMAT, &(sampleFormat));
  TIFFGetField(tiffds, TIFFTAG_BITSPERSAMPLE, &(dataSizeFileIn));
  TIFFGetField(tiffds, TIFFTAG_COMPRESSION, &(isCompression));

  tiffMetaData->isCompressed = isCompression;
  tiffMetaData->sampleFormat = sampleFormat;
  tiffMetaData->dataSizeFileIn = dataSizeFileIn/8;

  int tileRet = 1;
  int stripRet = 1;
  tileRet &= TIFFGetField(tiffds, TIFFTAG_TILEWIDTH, &(tiffMetaData->blocks_x_size));
  tileRet &= TIFFGetField(tiffds, TIFFTAG_TILELENGTH, &(tiffMetaData->blocks_y_size));
  tileRet &= TIFFGetField(tiffds, TIFFTAG_TILEOFFSETS, &offset);
  stripRet&=TIFFGetField(tiffds, TIFFTAG_ROWSPERSTRIP, &(tiffMetaData->blocks_y_size));
  stripRet&=TIFFGetField(tiffds, TIFFTAG_STRIPOFFSETS, &offset);

  if (tileRet == 1) {
    tiffMetaData->tiles_across = (tiffMetaData->x_size+ tiffMetaData->blocks_x_size - 1) / tiffMetaData->blocks_x_size;
    tiffMetaData->tiles_down = (tiffMetaData->y_size+ tiffMetaData->blocks_y_size - 1) / tiffMetaData->blocks_y_size;
    tiles_per_image = tiffMetaData->tiles_across * tiffMetaData->tiles_down;

    tiffMetaData->tile_offsets = new int64_t[tiles_per_image];
    if (tiffMetaData->tile_offsets == NULL) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error: unable to initialize memory for tile_offsets" <<endl;
      return false;
    }
    tileRet = TIFFGetField(tiffds, TIFFTAG_TILEOFFSETS, &tiff_offsets);
    memcpy(tiffMetaData->tile_offsets,tiff_offsets,sizeof(int64_t) * tiles_per_image);
  }

  if (stripRet==1) {
    tiffMetaData->tileOrStrip = STRIPTYPE;
    tiffMetaData->blocks_x_size = tiffMetaData->x_size;
    tiffMetaData->tiles_across = (tiffMetaData->x_size+ tiffMetaData->blocks_x_size - 1) / tiffMetaData->blocks_x_size;
    tiffMetaData->tiles_down = (tiffMetaData->y_size+ tiffMetaData->blocks_y_size - 1) / tiffMetaData->blocks_y_size;
    tiles_per_image = tiffMetaData->tiles_across * tiffMetaData->tiles_down;

    tiffMetaData->tile_offsets = new int64_t[tiles_per_image];
    if (tiffMetaData->tile_offsets == NULL) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error: unable to initialize memory for tile_offsets" <<endl;
      return false;
    }
    stripRet = TIFFGetField(tiffds, TIFFTAG_TILEOFFSETS, &tiff_offsets);
    memcpy(tiffMetaData->tile_offsets,tiff_offsets,sizeof(int64_t) * tiles_per_image);
  }

  if (stripRet != 1 && tileRet!=1) {
    tiffMetaData->tileOrStrip =STRIPTYPE;
    int ret;
    ret = TIFFGetField(tiffds, TIFFTAG_STRIPOFFSETS, &offset);
    if (ret != 1) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error occurred when reading strip offsets!" <<endl;
      return false;
    }
    tiffMetaData->first_strip_offset = *offset;
    tiffMetaData->tiles_across = 1;
    tiffMetaData->tiles_down = tiffMetaData->y_size;
    tiffMetaData->blocks_y_size = tiffMetaData->x_size;
  }
  tiffMetaData->first_strip_offset = *offset;
  TIFFClose(tiffds);
  return true;
}

bool gTiffMetaData::setOutputFile(char*fileName) {
  if(fileName==NULL){
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: output filename is null (" <<fileName<<")"<<endl;
    return false;
  }
  strcpy(outputGTiffFile, fileName);
  return true;
}

bool  gTiffMetaData::writeMetaData() {
 
  GDALAllRegister();
  GDALDriver *driver = GetGDALDriverManager()->GetDriverByName("GTiff");
  if (driver == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: unable to load the GTiff driver." <<endl;
    return false;
  }
  /* guan edit
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);*/
 
  string ts_x_str = intToString(tiffMetaData->blocks_x_size);
  string ts_y_str = intToString(tiffMetaData->blocks_y_size);
  char **options = NULL;
  options = CSLSetNameValue(options, "INTERLEAVE", "PIXEL");
  options = CSLSetNameValue(options, "BIGTIFF", "YES");
  options = CSLSetNameValue(options, "TILED", "YES");
  options = CSLSetNameValue(options, "COMPRESS", "NONE");
  options = CSLSetNameValue(options, "BLOCKXSIZE", ts_x_str.c_str());
  options = CSLSetNameValue(options, "BLOCKYSIZE", ts_y_str.c_str());
  options = CSLSetNameValue(options, "SPARSE_OK", "YES");
  GDALDataset *output =driver->Create(outputGTiffFile,tiffMetaData->x_size,tiffMetaData->y_size ,tiffMetaData->band_count,tiffMetaData->band_type,options);
  CSLDestroy(options);
  if (output == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: Failed to create the output raster file (" << outputGTiffFile << ")" \
        <<endl;
    return false;
  }
  GDALRasterBand* out_band = output->GetRasterBand(1);
  out_band->SetNoDataValue(tiffMetaData->no_data);
  output->SetGeoTransform(tiffMetaData->geoTransform);

  OGRSpatialReference out_sr;
  char *wkt = NULL;
  out_sr.SetFromUserInput(tiffMetaData->projectionRef.c_str());
  out_sr.exportToWkt(&wkt);

  output->SetProjection(wkt);
  OGRFree(wkt);
  GDALClose(output);

  MPI_File fh;
  int rc = MPI_File_open(MPI_COMM_SELF,
                         outputGTiffFile,
                         MPI_MODE_RDWR,
                         MPI_INFO_NULL,
                         &fh);
  MPI_File_set_atomicity(fh, 0);
  if (rc != MPI_SUCCESS) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: MPI_File failed to open file (" << outputGTiffFile << ")" <<endl;
    return false;
  }

  bool big_endian = false;   
  uint8_t endian_flag[2] = { 0x49, 0x49 };
  MPI_File_read_at(fh, 0, endian_flag, 2, MPI_BYTE, MPI_STATUS_IGNORE);
  if (endian_flag[0] == 0x4d) {
    big_endian = true;
  }

  uint8_t version_data[2];
  int16_t version = 0;
  MPI_File_read_at(fh, 2, version_data, 2, MPI_BYTE, MPI_STATUS_IGNORE);
  version = parse_int16(version_data, big_endian);

  if (version != 0x002B) {
    return false;
  }


  int64_t doffset = 0;
  doffset = read_int64(fh, 8, big_endian);

  int64_t entry_count = read_int64(fh, doffset, big_endian);
  int64_t entry_offset = doffset + sizeof(int64_t);

  int64_t tile_count = 0;
  tiffMetaData->band_type_size = GDALGetDataTypeSize(tiffMetaData->band_type)/8;

  const int64_t tile_size_bytes = tiffMetaData->blocks_x_size * tiffMetaData->blocks_y_size *
                                  tiffMetaData->band_count * tiffMetaData->band_type_size;
  int64_t first_tile_offset = 0;

  for (int64_t i = 0; i < entry_count; ++i) {
    uint8_t tag_buffer[2];
    MPI_File_read_at(fh,
                     entry_offset,
                     tag_buffer,
                     2,
                     MPI_BYTE,
                     MPI_STATUS_IGNORE);
    int16_t entry_tag = parse_int16(tag_buffer, big_endian);
    int64_t element_count = read_int64(fh, entry_offset+4, big_endian);
    int64_t entry_data = read_int64(fh, entry_offset+12, big_endian);

    if (entry_tag == TIFFTAG_TILEOFFSETS) {
      MPI_Offset first_offset;
      MPI_File_get_size(fh, &first_offset);
      first_tile_offset = first_offset;
      tile_count = element_count;

      for (int64_t j = 0; j < element_count; ++j) {
        write_int64(fh,entry_data+(sizeof(int64_t)*j),
                    first_offset+(tile_size_bytes*j),
                    big_endian);
      }
    }
    else if (entry_tag == TIFFTAG_TILEBYTECOUNTS) {
      for (int64_t j = 0; j < element_count; ++j) {
        write_int64(fh, entry_data+(sizeof(int64_t)*j),
                    tile_size_bytes, big_endian);
      }
    }
    entry_offset += 20;
  }
  uint8_t buffer = 0;
  int64_t file_size = (tile_count * tile_size_bytes) + first_tile_offset;
  MPI_File_write_at(fh, file_size-1, &buffer, 1, MPI_BYTE, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);

  return true;
}

string gTiffMetaData::intToString(int output_tile_size)
{
  string ts_str;
  stringstream ss;
  ss << output_tile_size;
  ss >> ts_str;
  return ts_str;
}



int64_t gTiffMetaData::read_int64(MPI_File&fh, int64_t offset, bool big_endian) {
  uint8_t read_buffer[8];
  MPI_File_read_at(fh,
                   offset,
                   read_buffer,
                   sizeof(int64_t),
                   MPI_BYTE,
                   MPI_STATUS_IGNORE);
  return parse_int64(read_buffer, big_endian);
}

int gTiffMetaData::write_int64(MPI_File&fh,
                               int64_t offset,
                               int64_t value,
                               bool big_endian) {
  uint8_t buffer[8];
  export_int64(value, buffer, big_endian);
  MPI_File_write_at(fh,
                    offset, buffer,
                    sizeof(int64_t),
                    MPI_BYTE,
                    MPI_STATUS_IGNORE);
  return 0;
}


int64_t gTiffMetaData::parse_int64(uint8_t *buffer, bool big_endian) {
  int64_t result = 0;
  int64_t temp = 0;

  if (big_endian) {
    temp = buffer[0];
    result |= temp<<56;
    temp = buffer[1];
    result |= temp<<48;
    temp = buffer[2];
    result |= temp<<40;
    temp = buffer[3];
    result |= temp<<32;
    temp = buffer[4];
    result |= temp<<24;
    temp = buffer[5];
    result |= temp<<16;
    temp = buffer[6];
    result |= temp<<8;
    temp = buffer[7];
    result |= temp<<0;
  }
  else {
    temp = buffer[7];
    result |= temp<<56;
    temp = buffer[6];
    result |= temp<<48;
    temp = buffer[5];
    result |= temp<<40;
    temp = buffer[4];
    result |= temp<<32;
    temp = buffer[3];
    result |= temp<<24;
    temp = buffer[2];
    result |= temp<<16;
    temp = buffer[1];
    result |= temp<<8;
    temp = buffer[0];
    result |= temp<<0;
  }

  return result;
}

int gTiffMetaData::export_int64(int64_t num, uint8_t *buffer, bool big_endian) {
  if (big_endian) {
    buffer[0] = (num>>56);
    buffer[1] = (num>>48);
    buffer[2] = (num>>40);
    buffer[3] = (num>>32);
    buffer[4] = (num>>24);
    buffer[5] = (num>>16);
    buffer[6] = (num>>8);
    buffer[7] = (num>>0);
  }
  else {
    buffer[7] = (num>>56);
    buffer[6] = (num>>48);
    buffer[5] = (num>>40);
    buffer[4] = (num>>32);
    buffer[3] = (num>>24);
    buffer[2] = (num>>16);
    buffer[1] = (num>>8);
    buffer[0] = (num>>0);
  }
  return 0;
}

bool gTiffMetaData::getMetaDataRaster()
{
  GDALAllRegister();
  /* guan Edit
  int _rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&_rank);
  */
  GDALDataset *ds = static_cast<GDALDataset*>(GDALOpen(outputGTiffFile, GA_ReadOnly));
  if (ds == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: unable to load the GTiff driver" << endl;
    return false;
  }

  tiffMetaData->x_size = ds->GetRasterXSize();
  tiffMetaData->y_size = ds->GetRasterYSize();
  tiffMetaData->band_count = ds->GetRasterCount();
  tiffMetaData->band_type = ds->GetRasterBand(1)->GetRasterDataType();
  tiffMetaData->band_type_size = GDALGetDataTypeSize(tiffMetaData->band_type)/8;
  tiffMetaData->first_strip_offset = -1;
  tiffMetaData->blocks_x_size = tiffMetaData->x_size;
  tiffMetaData->blocks_y_size = tiffMetaData->y_size;
  tiffMetaData->dataSizeObj = getDataObj(tiffMetaData->band_type);
  GDALClose(ds);

  TIFF *tiffds = TIFFOpen(outputGTiffFile, "r");
  if (tiffds == NULL) {
    cerr << __FILE__ << " " << __FUNCTION__ \
        << " Error: libtiff unable to open file (" << outputGTiffFile << ")" \
        << endl;
    return false;
  }
  int64_t *offset = NULL;
  int64_t tiles_per_image = -1;
  int64_t *tiff_offsets = NULL;

  int ret = 1;
  ret &= TIFFGetField(tiffds, TIFFTAG_TILEWIDTH, &(tiffMetaData->blocks_x_size));
  ret &= TIFFGetField(tiffds, TIFFTAG_TILELENGTH, &(tiffMetaData->blocks_y_size));
  ret &= TIFFGetField(tiffds, TIFFTAG_TILEOFFSETS, &offset);

  if (ret == 1) {
    tiffMetaData->tiles_across = (tiffMetaData->x_size + tiffMetaData->blocks_x_size - 1) /
                                 tiffMetaData->blocks_x_size;
    tiffMetaData->tiles_down = (tiffMetaData->y_size + tiffMetaData->blocks_y_size - 1) /
                               tiffMetaData->blocks_y_size;
    tiles_per_image = tiffMetaData->tiles_across * tiffMetaData->tiles_down;
    tiffMetaData->tile_offsets = new int64_t[tiles_per_image];
    if (tiffMetaData->tile_offsets == NULL) {
      return false;
    }
    ret = TIFFGetField(tiffds, TIFFTAG_TILEOFFSETS, &tiff_offsets);
    memcpy(tiffMetaData->tile_offsets,
           tiff_offsets,
           sizeof(int64_t) * tiles_per_image);
  }
  else if (ret != 1) {
    ret = TIFFGetField(tiffds, TIFFTAG_STRIPOFFSETS, &offset);
    if (ret != 1) {
      cerr << __FILE__ << " " << __FUNCTION__ \
          << " Error occurred when reading strip offsets" <<endl;
      return false;
    }
    tiffMetaData->first_strip_offset = *offset;
    tiffMetaData->tiles_across = 1;
    tiffMetaData->tiles_down = tiffMetaData->y_size;
    tiffMetaData->blocks_y_size = tiffMetaData->x_size;
  }
  tiffMetaData->first_strip_offset = *offset;
  TIFFClose(tiffds);
  return true;
}

int gTiffMetaData::getDataObj(GDALDataType band_type) {
  if (band_type == GDT_Byte){
    return sizeof (unsigned char);
  }
  else if (band_type == GDT_UInt16){
    return sizeof (unsigned short);
  }
  else if (band_type == GDT_UInt32){
    return sizeof (unsigned int);
  }
  else if (band_type == GDT_Int16){
    return sizeof (short);
  }
  else if (band_type == GDT_Int32){
    return sizeof (int);
  }
  else if (band_type== GDT_Float32){
    return sizeof (float);
  }else if(band_type == GDT_Float64){
	return sizeof (double);
  }
  return 1;
}

int16_t gTiffMetaData::parse_int16(uint8_t *buffer, bool big_endian) {
  if (big_endian)
    return (buffer[1]<<0 | buffer[0]<<8);
  else
    return (buffer[0]<<0 | buffer[1]<<8);
}
