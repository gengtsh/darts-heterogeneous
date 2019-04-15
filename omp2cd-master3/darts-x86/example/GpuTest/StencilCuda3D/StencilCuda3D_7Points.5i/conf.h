#define TILE_SIZE 1 //In elements 
#define MAX_BLOCK_DIM 8 
#define GRID_TILE_X 500
#define GRID_TILE_Y 20  


#define NUM_THREADS 256 
//TILE37_X*TILE37_X = NUM_THREADS
#define GRID_TILE37_X 16 
#define GRID_TILE37_Y 16
#define GRID_TILE37_Z 16//10 is only for test, should use #<NUM_THREADS 
#define GRID_TILECPU_X 16
#define GRID_TILECPU_Y 16
#define GRID_TILECPU_Z 10

#define DIM 3

#define NUMTILESMAX 20
#define NUMTILES 3
#define NUMASTREAMS 5
#define NUMSSTREAMS 4
#define NUMPAR 6
#define NUMPARMAX 50 
#define NUMSTREAMMAX NUMPARMAX*2
//numchunkalloc >= (NumTile*Tile+2)
#define NUMCHUNKALLOCX 80
#define NUMCHUNKALLOCY 80
#define NUMCHUNKALLOCZ 80

#define INCX 10
#define INCY 10
#define INCZ 10

#define NUMREC 10 //record number: 10

#define GPU_AVMEM_SIZE 2
//GPU available memory size equals to GPU_AVMEM_SIZE*GB

#define HALO 1  //TODO Not everything has this value in mind. Please update properly


