# pRPL
parallel Raster Processing Library (pRPL) is a C++ programming library that provides easy-to-use interfaces to parallelize raster/image processing algorithms.

<b>1. To Compile</b> <br>
Note: this version is not a final release, and some components are still under testsing. The program has been tested on  Scientific Linux 6.7 operation system, compiled using g++ 4.9, OpenMPI 1.9, GDAL 1.9, and LibTIFF 4.0. The makefile (i.e., make_pAspect) will compile a demonstration program, pAspect, which is able to calcuate aspect and slope from DEM data in parallel. <br>
(1) Before compiling, make sure MPI, GDAL, and LibTIFF libraries have been installed. <br>
(2) Open <b>make_pAspect</b> and modify the lines that specify the locations of libraries. <br>
(3) Type 'make -f make_pAspect depend'. <br>
(4) Type 'make -f make_pAspect' to compile. <br>
After successful compilation, an executable file named <b>pAspect</b> will be generated.

<b>2. To Run</b><br>
There is a <b>nd_dem.tif</b> file, which is the DEM data at 1.5km resolution of North Dakota in the US. Note that this data is only for testing if the pAspect program works as expected, not for demonstrating the performance. <br>

<b>2.1 Usage:</b><br>
mpirun -np \<num_proc\> pAspect \<workspace\> \<input-demFilename\> \<num-row-subspaces\> \<num-col-subspaces\> \<task-farming(1/0)\> \<io-option(0/1/2/3/4/5)\> \<with-writer(1/0)\>  <br>
<b>workspace</b>: the directory where the input file is located and the output files will be written. <br>
<b>input-demFilename</b>: the input file in the GeoTIFF format, usually the DEM data. <br>
<b>num-row-subspaces</b>: the number of sub-domains in the Y axis, for domain decomposition. If num-row-subspaces > 1 and num-col-subspaces = 1, the domain is decomposed as row-wise; if num-row-subspaces = 1 and num-col-subspaces > 1, the domain is decomposed as column-wise; if both > 1, the domain is decomposed as block-wise. <br>
<b>num-col-subspaces</b>: the number of sub-domains in the X axis, for domain decomposition. <br>
<b>task-farming</b>: load-balancing option, either 0 or 1. if 0, static load-balancing; if 1, task farming. <br>
<b>io-option</b>: I/O option, ranges within [0, 5]. Option 0: GDAL-based centralized reading, no writing; Option 1: GDAL-based parallel reading, no writing; Option 2: pGTIOL-based parallel reading, no writing; Option 3: GDAL-based centralized reading and writing; Option 4: GDAL-based parallel reading and pseudo parallel writing; Option 5: pGTIOL-based parallel reading and parallel writing. <br>
<b>with-writer</b>: an option that specify whether a writer process will be used. If 0, no writer; if 1, use a writer. <br>

<b>2.2 Example:</b><br>
mpirun -np 8 ./pAspect ./ nd_dem.tif 8 1 0 5 0 <br>

<b>3 To Cite in Publications:</b><br>
- Miao, J.; Guan, Q.*; Hu, S. 2017. pRPL + pGTIOL: The marriage of a parallel processing library and a parallel I/O library for big raster data. Environmental Modelling and Software, 96: 347-360
- Guan, Q.; Zeng, W.; Gong, J.; and Yun, S. 2014. pRPL 2.0: improving the parallel raster processing library. Transactions in GIS, 18(S1): 25-52
- Guan, Q. and Clarke, K. C. 2010. A general-purpose parallel raster processing programming library test application using a geographic cellular automata model. International Journal of Geographical Information Science, 24(5): 695-722





