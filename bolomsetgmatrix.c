#include <stdio.h>
#ifdef EFIT3365
#include "gmatrix3365.h"
#include "gpos3365.h"
#endif
#ifdef EFIT6565
#include "gmatrix6565.h"
#include "gpos6565.h"
#endif
extern float *gmatrix;
extern int *gpos;
extern int gmatrix_len;
extern int gpos_len;
bolom_set_gmatrix(shot)
int shot;
{
#ifdef EFIT6565
    if(shot >= gmatrix_7_firstshot){
        gmatrix = gmatrix_7;
        gpos = gpos_7;
        gmatrix_len = gmatrix_7_len;
        gpos_len = gpos_7_len;
    } else if(shot >= gmatrix_6_firstshot){
        gmatrix = gmatrix_6;
        gpos = gpos_6;
        gmatrix_len = gmatrix_6_len;
        gpos_len = gpos_6_len;
    } else if(shot >= gmatrix_5_firstshot){
        gmatrix = gmatrix_5;
        gpos = gpos_5;
        gmatrix_len = gmatrix_5_len;
        gpos_len = gpos_5_len;
    } else if(shot >= gmatrix_4_firstshot){
        gmatrix = gmatrix_4;
        gpos = gpos_4;
        gmatrix_len = gmatrix_4_len;
        gpos_len = gpos_4_len;
    } else if(shot >= gmatrix_3_firstshot){
        gmatrix = gmatrix_3;
        gpos = gpos_3;
        gmatrix_len = gmatrix_3_len;
        gpos_len = gpos_3_len;
    } else if(shot >= gmatrix_2_firstshot){
        gmatrix = gmatrix_2;
        gpos = gpos_2;
        gmatrix_len = gmatrix_2_len;
        gpos_len = gpos_2_len;
    } else if(shot >= gmatrix_1_firstshot){
        gmatrix = gmatrix_1;
        gpos = gpos_1;
        gmatrix_len = gmatrix_1_len;
        gpos_len = gpos_1_len;
    } else {
        fprintf(stderr,"No matrix for specified shot = %d\n",shot);
        return(0);
    }
#endif
#ifdef EFIT3365

// efit grid 33x65 is assumed to be deprecated. Matrices no longer created.

    if(shot >= gmatrix_5_firstshot){
        gmatrix = gmatrix_5;
        gpos = gpos_5;
        gmatrix_len = gmatrix_5_len;
        gpos_len = gpos_5_len;
    } else if(shot >= gmatrix_4_firstshot){
        gmatrix = gmatrix_4;
        gpos = gpos_4;
        gmatrix_len = gmatrix_4_len;
        gpos_len = gpos_4_len;
    } else if(shot >= gmatrix_3_firstshot){
        gmatrix = gmatrix_3;
        gpos = gpos_3;
        gmatrix_len = gmatrix_3_len;
        gpos_len = gpos_3_len;
    } else if(shot >= gmatrix_2_firstshot){
        gmatrix = gmatrix_2;
        gpos = gpos_2;
        gmatrix_len = gmatrix_2_len;
        gpos_len = gpos_2_len;
    } else if(shot >= gmatrix_1_firstshot){
        gmatrix = gmatrix_1;
        gpos = gpos_1;
        gmatrix_len = gmatrix_1_len;
        gpos_len = gpos_1_len;
    } else {
        fprintf(stderr,"No matrix for specified shot = %d\n",shot);
        return(0);
    }
#endif
		return(1);
}
