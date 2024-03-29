#include <stdio.h>
#ifdef EFIT3365
#include "bound3365.h"
#endif
#ifdef EFIT6565
#include "bound6565.h"
#endif
#include "bolom.h"
extern float *gmatrix;
extern int *gpos;
extern int gmatrix_len;


#define PROJ (gpos[k])
#define CELL (gpos[j])
bolomproj(shot,im,proj)
int shot;
image im;
float *proj;
{

    register int i,j,k,l;
    int chord;
    int num_of_projs;
    int cells;
    int xlen,ylen;
    int len;

    bolom_set_gmatrix(shot);

    cells = XLEN * YLEN;
    num_of_projs = CHANS;
    len = gmatrix_len;

    memset(proj,0,num_of_projs * sizeof(float));
    for(i = 0,j=0,k=1;i < len;++i,j+=2,k+=2){
        if(CELL < cells && PROJ < num_of_projs && !bound_flag[CELL]){
            proj[PROJ] += gmatrix[i] * im[CELL];
        }
    }
}
bolombackproj(shot,im,proj)
int shot;
image im;
float *proj;
{

    register int i,j,k,l;
    int chord;
    int num_of_projs;
    int cells;
    int xlen,ylen;
    int len;

    bolom_set_gmatrix(shot);

    cells = XLEN * YLEN;
    num_of_projs = CHANS;
    len = gmatrix_len;

    memset(im,0,cells * sizeof(float));
    for(i = 0,j=0,k=1;i < len;++i,j+=2,k+=2){
        if(CELL < cells && PROJ < num_of_projs && !bound_flag[CELL]){
	    im[CELL] += gmatrix[i] * proj[PROJ];
        }
    }
}
bolomclip(shot,im,chans,nchans)
int shot;
image im;
int *chans,nchans;
{

    register int i,j,k,l;
    int chord;
    int num_of_projs;
    int cells;
    int xlen,ylen;
    int len;
		int clip[XLEN*YLEN];
		int chansused[CHANS];

    bolom_set_gmatrix(shot);

    cells = XLEN * YLEN;
    num_of_projs = CHANS;
    len = gmatrix_len;

    memset(clip,0,XLEN*YLEN*sizeof(int));
    memset(chansused,0,CHANS*sizeof(int));

/* Leonard chose to select region based on first 48 channels  */
		/*for(i = 0; i < nchans; ++i) chansused[chans[i]] = 1;*/
		for(i = 0; i < 48; ++i) chansused[i] = 1;

    for(i = 0,j=0,k=1;i < len;++i,j+=2,k+=2){
        if(CELL < cells && PROJ < CHANS && chansused[PROJ]){
            clip[CELL] = 1;
        }
    }

		for(i = 0; i < cells; ++i){
				if(clip[i] == 0){
					 im[i] = 0;
        }
		}
}
write_projdat(name,proj)
char *name;
float *proj;
{

    FILE *fefit;
    int i;
    fefit = fopen(name,"w");
    for(i = 0; i < CHANS; ++i){
        fprintf(fefit,"%d %g\n",i,proj[i]);
    }
    fclose(fefit);


}
write_imagedat(name,im)
char *name;
image im;
{

    FILE *fefit;
    int i;
    fefit = fopen(name,"w");
    fwrite(im,1,XLEN * YLEN * sizeof(float),fefit);
    fclose(fefit);


}
dump_imagedat(im)
image im;
{

    int i;
    for(i = 0; i < XLEN*YLEN; ++i)printf("%g\n",im[i]);


}
