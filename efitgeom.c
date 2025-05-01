#include "bolom.h"
#include <math.h>


void efitgeom(float *efit,float rax,float zax,float rxpt1,float zxpt1,float rxpt2,float zxpt2,float *psi_norm)
{
    int i,j;
    int ai,aj;
    float ig;
    float xinc,yinc;
    float slope,intercept;



    xinc = (XMAX - XMIN) / (float)(XLEN-1);
    yinc = (YMAX - YMIN) / (float)(YLEN-1);
    ai = (rax - XMIN) / xinc;
    aj = (zax - YMIN) / yinc;
    if(rxpt1 < XMIN || rxpt1 > XMAX ||
        zxpt1 < YMIN || zxpt1 > YMAX){
        rxpt1 = rax;
        zxpt1 = YMIN;
    }
    if(rxpt2 < XMIN || rxpt2 > XMAX ||
        zxpt2 < YMIN || zxpt2 > YMAX){
        rxpt2 = rax;
        zxpt2 = YMAX;
    }


    if(rxpt1 == rax) slope = 0.0;
    else{
        slope = (zax - zxpt1) / yinc / (rax - rxpt1) * xinc;
        intercept = (zxpt1 - (float)YMIN) / yinc - slope * (rxpt1 - (float)XMIN) / xinc;
    }
    for(j=0; j <= aj; ++j){
        if(slope == 0.0)ig = ai;
        else ig = ((float)j - intercept) / slope;
        for(i = 0; i < XLEN ; ++i){
            if((float)i < ig) psi_norm[i+j*XLEN] = -efit[i+j*XLEN];
            else psi_norm[i+j*XLEN] = efit[i+j*XLEN];
        }
    }

    if(rxpt2 == rax) slope = 0.0;
    else{ 
        slope = (zax - zxpt2) / yinc / (rax - rxpt2) * xinc;
        intercept = (zxpt2 - (float)YMIN) / yinc - slope * (rxpt2 - (float)XMIN) / xinc;
    }
    for(j = aj; j < YLEN; ++j){
        if(slope == 0.0)ig = ai;
        else ig = ((float)j - intercept) / slope;
        for(i = 0; i < XLEN ; ++i){
            if((float)i < ig) psi_norm[i+j*XLEN] = -efit[i+j*XLEN];
            else psi_norm[i+j*XLEN] = efit[i+j*XLEN];
        }
    }
}
