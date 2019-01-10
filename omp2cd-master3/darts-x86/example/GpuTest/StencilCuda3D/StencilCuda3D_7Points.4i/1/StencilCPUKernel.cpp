
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
//#include <math.h>
#include <cmath>
#include <cstring>
#include "conf.h"
#include "stencil.h"
#include "StencilCPUKernel.h"
//#define DEBUG_COPY

void computeInner_stencil2d_v2(double *dst,double *src,size_t n_rows,size_t n_cols){

	typedef double (*Array2D)[n_cols];
	Array2D DST = (Array2D) dst,
			SRC = (Array2D) src;

	for (size_t i = 1; i < n_rows; ++i) {
	    for (size_t j = 1; j < n_cols-1; ++j) {
	       // DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
	        DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 5.5;
		}
	}
}


void computeInner_stencil37(double *dst,double *src,size_t n_rows,size_t n_cols,size_t n_slices){

	typedef double (*Array3D)[n_rows][n_cols];
	Array3D DST = (Array3D) dst,
			SRC = (Array3D) src;

    for (size_t k = 1; k< n_slices;++k){
    	for (size_t i = 1; i < n_rows-1; ++i) {
            for (size_t j = 1; j < n_cols-1; ++j) {
                DST[k][i][j] = (SRC[k][i-1][j] + SRC[k][i+1][j] + SRC[k][i][j-1] + SRC[k][i][j+1] + SRC[k][i][j]+SRC[k-1][i][j]+SRC[k+1][i][j])/7.5;
            }
        }
    }
}


void computeBlock_stencil37(double *dst,double *src,size_t n_rows,size_t n_cols,size_t n_slices,size_t n_rows_ck,size_t n_cols_ck,size_t n_slices_ck){

	typedef double (*Array3D)[n_rows][n_cols];
	Array3D DST = (Array3D) dst,
			SRC = (Array3D) src;

    for (size_t k = 1; k< n_slices_ck;++k){
    	for (size_t i = 1; i < n_rows_ck; ++i) {
            for (size_t j = 1; j < n_cols_ck; ++j) {
                DST[k][i][j] = (SRC[k][i-1][j] + SRC[k][i+1][j] + SRC[k][i][j-1] + SRC[k][i][j+1] + SRC[k][i][j]+SRC[k-1][i][j]+SRC[k+1][i][j])/7.5;
            }
        }
    }
}


void copyColumns_stencil2d(double *dst,double *src,size_t n_rows,size_t n_cols,size_t stride){
	int n_cols_cp = std::ceil(1.0*n_cols/stride);

	typedef double (*Array2D)[n_cols];
	
	Array2D SRC = (Array2D) src;
	
	for (size_t i = 0;i<n_rows;++i){
		dst[i] = src[i*stride];
#ifdef DEBUG_COPY
		std::cout<<"shared_cols["<<i<<"] = "<< dst[i]<<std::endl;
#endif
	}
}


void copyRows_stencil2d(double *dst,double *src,size_t n_rows,size_t n_cols,size_t stride){
	int cpRow;
	int sz = std::ceil(1.0*n_rows/stride);
	for (int i = 0;i<sz;++i){
		cpRow = i*n_cols*stride;
		std::copy(src+cpRow, src+n_cols*2+cpRow, dst+i*n_cols*2);
	}

#ifdef DEBUG_COPY
		std::cout<<"shared_Rows:"<<std::endl;
	for(int i=0;i<sz;++i){	
		for(int j =0;j<2*n_cols;++j){
			std::cout<<dst[2*n_cols*i+j]<<",";
		}
		std::cout<<"\n";
	}
#endif

}
void calcNextEEL(int * arrncgEdge, int *arrncgEdgeLeft,int *arrnEdge, int *varEdge,int *arrncgEdgeMin,int id, int homeCnt, int visitCnt){
        //int newEdge;
        //int newEdgeLeft;
        //if(homeCnt>visitCnt){
        //    newEdge =  (arrncgEdge[id]+varEdge[id]);
        //}else{
        //    newEdge = (varEdge[id]>arrncgEdge[id])? arrncgEdge[id] : (arrncgEdge[id]-varEdge[id]);

        //}
        //newEdge = (newEdge > arrncgEdgeMin[id])? newEdge :arrncgEdgeMin[id] ;
        //arrncgEdge[id] = (newEdge>arrncgEdgeLeft[id])?(arrncgEdgeLeft[id]):(newEdge) ;
        //arrncgEdgeLeft[id] = arrncgEdgeLeft[id]-arrncgEdge[id];

        int newEdge;
        int newEdgeLeft;
        if(homeCnt>visitCnt){
            newEdge =  (arrncgEdge[id]+varEdge[id]);
        }else{
            newEdge = (varEdge[id]>arrncgEdge[id])? arrncgEdge[id] : (arrncgEdge[id]-varEdge[id]);

        }
        newEdge = (newEdge > arrncgEdgeMin[id])? newEdge :arrncgEdgeMin[id] ;
        newEdgeLeft = arrncgEdgeLeft[id] - newEdge;
        newEdge = (newEdgeLeft>arrncgEdgeMin[id])? newEdge : arrncgEdgeLeft[id] ;
        arrncgEdge[id] = (newEdge>arrncgEdgeLeft[id])?(arrncgEdgeLeft[id]):(newEdge) ;
        arrncgEdgeLeft[id] = arrncgEdgeLeft[id]-arrncgEdge[id];

}


void calcNextEP(int *arrncEdge, int *arrngEdge,int *arrncEdgeLeft, int *arrngEdgeLeft,int *cpos, int *gpos,int *arrnEdge, int *cvarEdge, int *gvarEdge,int *arrncEdgeMin,int *arrngEdgeMin,int cCnt,int gCnt, int nId,char *home){
/*
 *R(nId) level 3, R(nId+1) level 2,R(nId-1) level 1,
 *stability  3>2>1,  so, 1 change, then 2 change, then 3 change
 *homeCnt is current device, visitCnt is other device, e.g. current calculate cpu edge, homeCnt=cpucnt, visitCnt=gpucnt
 */
    int rnId = R(nId);
    int rnIdP1 = R(nId+1);
    int rnIdM1 = R(nId-1);
    
    int *arrnHomeEdge;
    int *arrnVisitEdge;
    int *arrnHomeEdgeLeft;
    int *arrnVisitEdgeLeft;
    int *arrnHomeEdgeMin;
    int *arrnVisitEdgeMin;
    int *homePos;
    int *visitPos;
    int *varEdge;
    int homeCnt;
    int visitCnt;
    int sign;
    
    if (strcmp(home,"cpu")==0){
        arrnHomeEdge  = arrncEdge;
        arrnVisitEdge = arrngEdge;
        arrnHomeEdgeLeft  = arrncEdgeLeft;
        arrnVisitEdgeLeft = arrngEdgeLeft;
        arrnHomeEdgeMin   = arrncEdgeMin;
        arrnVisitEdgeMin  = arrngEdgeMin;
        homePos  = cpos;
        visitPos = gpos;
        varEdge = cvarEdge;
        homeCnt = cCnt;
        visitCnt = gCnt;
        sign = 1;
    }else if (strcmp(home,"gpu")==0){
        arrnHomeEdge  = arrngEdge;
        arrnVisitEdge = arrncEdge;
        arrnHomeEdgeLeft  = arrngEdgeLeft;
        arrnVisitEdgeLeft = arrncEdgeLeft;
        arrnHomeEdgeMin   = arrngEdgeMin;
        arrnVisitEdgeMin  = arrncEdgeMin;
        homePos  = gpos;
        visitPos = cpos;
        varEdge = gvarEdge;
        homeCnt = gCnt;
        visitCnt = cCnt;
        sign = -1;
    }

    int newEdge;

    if(arrnHomeEdgeLeft[rnIdM1]>0){// check whether HomeEdge[R(nId-1)] finish or not
        if(sign ==1){ //cpu position
            homePos[rnIdM1] = homePos[rnIdM1] + arrnHomeEdge[rnIdM1];
        }
        calcNextEEL(arrnHomeEdge, arrnHomeEdgeLeft,arrnEdge, varEdge,arrnHomeEdgeMin, rnIdM1, homeCnt, visitCnt);
        if(sign ==-1){ //gpu position
            homePos[rnIdM1] = homePos[rnIdM1] - arrnHomeEdge[rnIdM1];
        }
        if((homePos[rnId] == visitPos[rnId]) &&(homePos[rnIdP1]==visitPos[rnIdP1])){//check whether home and visit in the same rnId, rnIdP1
            arrnVisitEdgeLeft[rnIdM1] = arrnHomeEdgeLeft[rnIdM1];
        }
    }else{//arrnHomeEdgeLeft[R(nId-1)]==0
        
        if(arrnHomeEdgeLeft[rnIdP1]>0){// check whether HomeEdge[R(nId+1)] finish or not
            if(sign ==1){//cpu position
                homePos[rnIdP1] = homePos[rnIdP1] + arrnHomeEdge[rnIdP1];
                homePos[rnIdM1] = 0 ;
            }
            
            arrnHomeEdgeLeft[rnIdM1] = arrnEdge[rnIdM1] - arrnHomeEdge[rnIdM1];
            calcNextEEL(arrnHomeEdge, arrnHomeEdgeLeft,arrnEdge, varEdge,arrnHomeEdgeMin, rnIdP1, homeCnt, visitCnt);
            
            if(sign == -1){//gpu position
                homePos[rnIdP1] = homePos[rnIdP1]  - arrnHomeEdge[rnIdP1];
                homePos[rnIdM1] = arrnEdge[rnIdM1] - arrnHomeEdge[rnIdM1] ;
            }
            
            if(homePos[rnId]==visitPos[rnId]){ //home and visit in the same rnId
               arrnVisitEdgeLeft[rnIdP1] = arrnHomeEdgeLeft[rnIdP1]; 
            }

        }else{//arrnHomeEdgeLeft[R(nId+1)]==0
            
            if(arrnHomeEdgeLeft[rnId]>0){// check whether Edge[R(nId)] finish or not
                if(sign == 1){//cpu position
                    homePos[rnId] = homePos[rnId] + arrnHomeEdge[rnId];
                    homePos[rnIdP1] = 0 ;
                    homePos[rnIdM1] = 0 ;
                }
                arrnHomeEdgeLeft[rnIdM1] = arrnEdge[rnIdM1] - arrnHomeEdge[rnIdM1];
                arrnHomeEdgeLeft[rnIdP1] = arrnEdge[rnIdP1] - arrnHomeEdge[rnIdP1];
                calcNextEEL(arrnHomeEdge, arrnHomeEdgeLeft,arrnEdge, varEdge,arrnHomeEdgeMin, rnId, homeCnt, visitCnt);
                if(sign == -1){//gpu position
                    homePos[rnId] = homePos[rnId] - arrnHomeEdge[rnId];
                    homePos[rnIdP1] = arrnEdge[rnIdP1] - arrnHomeEdge[rnIdP1] ;
                    homePos[rnIdM1] = arrnEdge[rnIdM1] - arrnHomeEdge[rnIdM1] ;
                }

            }else{//arrnHomeEdgeLeft[R(nId)]==0
                
                if(arrnVisitEdgeLeft[rnIdP1]>0){// check whether visit Edge[rnIdP1] finish or not
                    //home and visit are in same rnId 
                    arrnHomeEdgeLeft[rnIdM1] = arrnEdge[rnIdM1] - arrnHomeEdge[rnIdM1];
                    arrnHomeEdgeLeft[rnIdP1] = arrnVisitEdgeLeft[rnIdP1]; 
                    homePos[rnId] = visitPos[rnId];
                    homePos[rnIdM1] = (sign == 1)? 0 : arrnEdge[rnIdM1] - arrnHomeEdge[rnIdM1] ;
                    arrnHomeEdge[rnId] = arrnVisitEdge[rnId];
                    calcNextEEL(arrnHomeEdge, arrnHomeEdgeLeft,arrnEdge, varEdge,arrnHomeEdgeMin, rnIdP1, homeCnt, visitCnt);
                    homePos[rnIdP1] = (sign == 1)? 0 : arrnEdge[rnIdP1] - arrnHomeEdge[rnIdP1] ;
                    arrnVisitEdgeLeft[rnIdP1] = arrnHomeEdgeLeft[rnIdP1]; 
                
                }else{ //arrnVisitEdge[rnIdP1] finish
                    
                    if(arrnVisitEdgeLeft[rnIdM1]>0){ // check whether visit Edge[rnIdM1] finish or not, 
                        //home and visit are in same rnId,rnIdP1
                        arrnHomeEdgeLeft[rnIdM1] = arrnVisitEdgeLeft[rnIdM1];
                        arrnHomeEdge[rnIdP1] = arrnVisitEdge[rnIdP1];
                        homePos[rnId] = visitPos[rnId];
                        homePos[rnIdP1] = visitPos[rnIdP1]; 
                        calcNextEEL(arrnHomeEdge, arrnHomeEdgeLeft,arrnEdge, varEdge,arrnHomeEdge, rnIdM1, homeCnt, visitCnt);
                        homePos[rnIdM1] = (sign == 1)? 0 : arrnEdge[rnIdM1] - arrnHomeEdge[rnIdM1] ;
                        arrnVisitEdgeLeft[rnIdM1] = arrnHomeEdgeLeft[rnIdM1]; 
                    }
                }
            }
        }
    }
    //both cpu and gpu have the same arrnEdgeLeft[R(nId)]
    arrncEdgeLeft[rnId] = arrnHomeEdgeLeft[rnId];
    arrngEdgeLeft[rnId] = arrnHomeEdgeLeft[rnId];
}


