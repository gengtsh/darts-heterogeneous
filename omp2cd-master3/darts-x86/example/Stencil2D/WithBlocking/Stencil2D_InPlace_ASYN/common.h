#ifndef STENCIL_COMMON
#define STENCIL_COMMON

//#define VERIFICATION 
//#define VERIFICATION_PRINT 
//#define STENCIL2DDEBUG

#include<cmath>

#define DIM 2

#define GRID_TILECPU_X 48 
#define GRID_TILECPU_Y 48

#define kb(x) (size_t (x)<<10)
#define mb(x) (size_t (x)<<20)
#define gb(x) (size_t (x)<<30)
#define KB kb(1)
#define MB mb(1)
#define XMB mb(20)
#define GB gb(1)

#define SWAP(type,left,right) do { \
    type tmp = left;               \
    left     = right;              \
    right    = tmp;                \
} while(0)

inline void swap_ptr(void** left, void** right) {

	void* tmp = *left;
    *left     = *right;
    *right    = tmp;
}

#define SWAP_PTR(left,right) swap_ptr((void**)left,(void**)right)



struct SRECORD{
   int chunk[DIM];
   size_t tStart;
   size_t tEnd;
   size_t tExe;
};


inline int fRR(int id,int dim){
    if (id>=dim) {
        return id-dim;
    }else if(id<0){
        return dim+id;
    }else{
        return id;
    }
}

#define R(ID) fRR((int)ID,(int)DIM)

inline int fRR2(int id, int dim){
    return id%dim;
}

#define RT(ID) fRR2((int)ID,(int)NUMREC)

inline void chooseSmaller(int *arrnSmaller,int *arrnFirst, int *arrnSecond, int offset1,int offset2,int offset3,int offset4){
   
    for(int i=0; i< DIM; ++i){
        arrnSmaller[i]=((arrnFirst[i]+offset1)>(arrnSecond[i]+offset2))? (arrnSecond[i]+offset3):(arrnFirst[i]+offset4);
    }
}

inline void cmpsAndassign(int *arrnNew,int *arrnFirst, int *arrnSecond, int *arrnThird, int*arrnFourth, int offset1,int offset2,int offset3,int offset4){
   
    for(int i=0; i< DIM; ++i){
        arrnNew[i]=((arrnFirst[i]+offset1)<(arrnSecond[i]+offset2))? (arrnThird[i]+offset3):(arrnFourth[i]+offset4);
    }
}


inline void calcarrnDivCeil(int *arrnDivCeil, int *arrnFirst,int *arrnSecond){
    for(int i=0;i<DIM;++i){
        arrnDivCeil[i] = std::ceil(1.0*arrnFirst[i]/arrnSecond[i]);
    }
}

inline void calcarrnCpuEdgeLeft(int *arrnEdgeLeft, int *arrlEdge, int *arrnPos, int *arrnEdge ){
    for(int i=0;i<DIM;++i){
        arrnEdgeLeft[i] = arrlEdge[i] - arrnPos[i] - arrnEdge[i];
    }
}

inline void calcarrnGpuEdgeLeft(int *arrnEdgeLeft, int *arrlEdge, int *arrnPos, int *arrnEdge ){
    for(int i=0;i<DIM;++i){
        arrnEdgeLeft[i] = arrnPos[i] ;
    }
}


inline void calcEdgeChunk(int *arrnEdgeChunk,int *pos,int *arrnChunk, int *arrnEdge, int*arrnChunk2,int offset1,int offset2,int offset3,int offset4){
   
    for(int i=0; i< DIM; ++i){
        arrnEdgeChunk[i]=(((pos[i]+offset1)*arrnChunk[i])>=(arrnEdge[i]+offset2))? (arrnEdge[i]-pos[i]*arrnChunk[i]+offset3):(arrnChunk2[i]+offset4);
    }
}

inline void calcEdge(int *arrnEdgeNew,int *pos,int *arrnEdge1, int *arrnEdge0,int offset1,int offset2,int offset3,int offset4){
   
    for(int i=0; i< DIM; ++i){
        arrnEdgeNew[i]=((pos[i]+arrnEdge1[i]+offset1)>=(arrnEdge0[i]+offset2))? (arrnEdge0[i]-pos[i]):(arrnEdge1[i] +offset4);
    }
}


inline void setarrnValue(int *arrn,int value){
    for(int i=0; i< DIM; ++i){
        arrn[i] = value;
    }
}

inline void setarrn1Fromarrn2(int *arrn1,int *arrn2){
    for(int i=0; i< DIM; ++i){
        arrn1[i] = arrn2[i];
    }
}

inline void calcarrn1Fromarrn2Marrn3(int *arrn1,int *arrn2,int *arrn3){
    for(int i=0; i< DIM; ++i){
        arrn1[i] = arrn2[i]*arrn3[i];
    }
}

inline void calcarrn1Fromarrn2Minusarrn3(int *arrn1,int *arrn2,int *arrn3){
    for(int i=0; i< DIM; ++i){
        arrn1[i] = arrn2[i]-arrn3[i];
    }
}

inline bool checkarrnEqValue(int *arrn,int value){
    for(int i=0; i< DIM; ++i){
        if(arrn[i] != value){
            return false;
        }
    }
    return true;
}


#endif // STENCIL_COMMON
