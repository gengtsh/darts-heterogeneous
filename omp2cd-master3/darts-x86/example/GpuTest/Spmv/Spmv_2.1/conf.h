#ifndef SPMV_CONF_H_
#define SPMV_CONF_H_


//#define DARTS_DEBUG
//#define DARTS_DEBUG_CE
#define DARTS_RECORD

//#define DARTS_DEBUG_VAL

//#define CUDA_RECORD

#define CUDA_V10  0
#define CUDA_V9   1

static const int BLOCK_SIZE = 128;
static const int WARP_SIZE = 32;
static const int NSTREAM   = 128;

#endif
