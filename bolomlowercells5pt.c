
/* 20110323 whm changed lsei to define, defaults to lsei_  */
/* 20031031 tbt changed lsei to lsei_  */
#if defined(nounderscore)
#define LSEI lsei
#else
#define LSEI lsei_
#endif


#define PREC float
#define WS_LEN 256000
#define IP_LEN 8192
#include <stdio.h>
/*#include <VIEW/user.h>*/
#include <nrutil.h>
#include <math.h>
#define EMIN 0.0
#define EMAX 20.0
static int debug  = 0;
#include "bolom.h"

#include <math.h>
#include <stdlib.h>
#include <nrutil.h>
#include <signal.h>
#include <stdio.h>
#include "boundext.h"
#include <setjmp.h>
static jmp_buf sjbuf ;
static int ilxp;
static float xinc,yinc;
float *gmatrix;
int *gpos;
int gmatrix_len;
int gpos_len;
#define PROJ (gpos[k])
#define CELL (gpos[j])
int grad[XLEN*YLEN];
image psi_norm;

static void (*prev_handler)();
static float *undump;
static float *insigma;
static float *fitpts;
static int nfitpts;
static image inefit;
static int done;
static int iuix,ouix;


bolomfit_lower_cells_5pt(shot,chans,nchans,inproj,sigma,zinner,zouter,minprivate,mincore,maxinner,maxouter,drweight,dzweight,efit,rax,zax,
rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,fret)
unsigned int shot;
float *inproj,*sigma,*fitproj,*fret;
float zinner,zouter,minprivate,mincore,maxinner,maxouter;
float drweight,dzweight;
image efit,fitimage;
int *chans,nchans;
float rax,zax,rxpt1,zxpt1,rxpt2,zxpt2;

{

    float *p,**xi;
    int iter;
    int splcore();
    FILE *f;
    int i,j,k;
    PREC c1,c2,c3,c4;
    PREC *amat,*cmat,*imat,*timat;
    PREC sum;
    PREC *work;
    int na,nma,eqn,toteqn;
    int l_m,l_n,l_p,l_lda,l_ldb;
    int lwork,info;
    float basisel(),basisv();
    int lda,numunk,numlow;
    int iparm;
    float *g,*h,*w;
    int mdw,me,ma,mg;
    int *ip;
    float *ws;
    float prgopt[100],rnorme,rnorml;
    int mode;
    int allbound;
    int here;


    PREC *u;
    PREC *tz,*uz;
    PREC *sol;

    struct sigvec vec,ovec;
    efitgeom(efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,psi_norm);
    if(!bolom_set_gmatrix(shot))return(0);


    done = 0;
    setjmp(sjbuf);
    if(done)return;
    done = 1;

    vec.sv_handler = splcore;
    vec.sv_mask = 0xffff;
    vec.sv_onstack = 0;
    vec.sv_flags = 0;
    sigvector(SIGINT,&vec,&ovec);
    prev_handler = ovec.sv_handler;
    sleep(10);

    /*write_projdat("coreorig.dat",inproj);*/

    xinc = (XMAX - XMIN) / (float)(XLEN-1);
    yinc = (YMAX - YMIN) / (float)(YLEN-1);

    if(rxpt1 < XMIN || rxpt1 > XMAX ||
        zxpt1 < YMIN || zxpt1 > YMAX){
        ilxp = 0;
    } else ilxp = (zxpt1 - YMIN) / yinc;
    ilxp *= XLEN;

    if(zinner < YMIN || zinner > YMAX){
        iuix = 0;
    } else iuix = (zinner - YMIN) / yinc;
    iuix *= XLEN;

    if(zouter < YMIN || zouter > YMAX){
        ouix = 0;
    } else ouix = (zouter - YMIN) / yinc;
    ouix *= XLEN;

    memset(grad,0L,sizeof(grad));
    here = 0;


    printf("ilxp = %d  iuix = %d  ouix = %d\n",ilxp,iuix,ouix);

    sigvector(SIGINT,&ovec,0);





    na = XLEN * YLEN;
    /*printf("na:%d nchans:%d\n",na,nchans);*/
    numlow = 0;

/*  Do a flat field back projection so we exclude cells with no measured response */
    for(i = 0; i < CHANS; ++i) fitproj[i] = 1;
    bolombackproj(shot,fitimage,fitproj);

    for(i = 0;i < ilxp;++i){
        if(bound_flag[i] || fitimage[i] == 0.0)continue;
        if(psi_norm[i] <= 0 && efit[i] >= minprivate &&
            efit[i] <= maxinner){
            grad[i] = numlow + 1.0;
            ++numlow;
        }
        if(psi_norm[i] > 0 && efit[i] >= minprivate &&
            efit[i] <= maxouter){
            grad[i] = numlow + 1.0;
            ++numlow;
        }
    }
    for(i = ilxp;i < iuix;++i){
        if(bound_flag[i] || fitimage[i] == 0.0)continue;
        if(psi_norm[i] <= 0 && efit[i] >= mincore &&
            efit[i] <= maxinner){
            grad[i] = numlow + 1.0;
            ++numlow;
        }
    }
    for(i = ilxp;i < ouix;++i){
        if(bound_flag[i] || fitimage[i] == 0.0)continue;
        if(psi_norm[i] > 0 && efit[i] >= mincore &&
            efit[i] <= maxouter){
            grad[i] = numlow + 1.0;
            ++numlow;
        }
    }
    lda = numlow;
    toteqn = 2*lda + nchans;
    printf("numlow = %d lda = %d toteqn = %d\n",numlow,lda,toteqn);
    cmat = calloc(sizeof(PREC),na);

    imat = calloc(sizeof(PREC) * 2,lda * toteqn);
    timat = calloc(sizeof(PREC) * 2,lda * toteqn);

    amat = calloc(sizeof(PREC) * 2,na * (CHANS + 1));
    for(i = 0,j=0,k=1;i < gmatrix_len;++i,j+=2,k+=2){
        if(CELL < na && PROJ < CHANS ){
            amat[PROJ * na + CELL] = gmatrix[i];
        }
    }
    for(i=0;i<nchans;++i){
        for(j = 0,k = 0;k < na;++k){
            if(grad[k]){
                imat[i * lda + j] = amat[chans[i]*na+k];
                ++j;
            }
        }
    }
    for(j=0,k=0;k < na;++k){
        if(grad[k] > 0.0){
            imat[(j+nchans)*lda+j] += -4 * drweight;
            if((k%XLEN)>0 && grad[k-1] > 0.0){
                imat[(j+nchans)*lda+(j-1)] += 1 * drweight;
            }
            if((k%XLEN)>1 && grad[k-2] > 0.0){
                imat[(j+nchans)*lda+(j-2)] += 0.5 * drweight;
            }
            if((k%XLEN)<(XLEN-1) && grad[k+1] > 0.0){
                imat[(j+nchans)*lda+(j+1)] += 1 * drweight;
            }
            if((k%XLEN)<(XLEN-2) && grad[k+2] > 0.0){
                imat[(j+nchans)*lda+(j+2)] += 0.5 * drweight;
            }
            ++j;
        }
    }
    for(j=0,k=0;k < na;++k){
        if(grad[k] > 0.0){
            imat[(j+nchans+numlow)*lda+j] += -4 * dzweight;
            if(k>XLEN && grad[k-XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[k-XLEN]-1] += 1 * dzweight;
            }
            if(k>(2*XLEN) && grad[k-2*XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[k-2*XLEN]-1] += 0.5 * dzweight;
            }
            if(k<(na-XLEN) && grad[k+XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[k+XLEN]-1] += 1 * dzweight;
            }
            if(k<(na-2*XLEN) && grad[k+2*XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[k+2*XLEN]-1] += 0.5 * dzweight;
            }
            ++j;
        }
    }
    for(j=0,k=0;k < na;++k){
        if(grad[k] > 0.0){
            imat[(j+nchans+numlow)*lda+j] += -8 * (dzweight + drweight);
            i = k - XLEN - 1;
            if(k>XLEN && (k%XLEN)>0 && grad[i] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[i]-1] += 
                    (dzweight + drweight) / sqrt(2.0);
            }
            i = k - XLEN + 1;
            if(k>XLEN && (k%XLEN)<(XLEN-1) && grad[i] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[i]-1] += 
                    (dzweight + drweight) / sqrt(2.0);
            }
            i = k + XLEN - 1;
            if(k<(na-XLEN) && (k%XLEN)>0 && grad[k+XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[i]-1] += 
                    (dzweight + drweight) / sqrt(2.0);
            }
            i = k + XLEN + 1;
            if(k<(na-XLEN) && (k%XLEN)<(XLEN-1) && grad[k+XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[i]-1] += 
                    (dzweight + drweight) / sqrt(2.0);
            }
            i = k - 2*XLEN - 2;
            if(k>2*XLEN && (k%XLEN)>1 && grad[i] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[i]-1] += 
                    (dzweight + drweight) / sqrt(8.0);
            }
            i = k - 2*XLEN + 2;
            if(k>2*XLEN && (k%XLEN)<(XLEN-2) && grad[i] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[i]-1] += 
                    (dzweight + drweight) / sqrt(8.0);
            }
            i = k + 2*XLEN - 2;
            if(k<(na-2*XLEN) && (k%XLEN)>1 && grad[k+XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[i]-1] += 
                    (dzweight + drweight) / sqrt(8.0);
            }
            i = k + XLEN + 1;
            if(k<(na-2*XLEN) && (k%XLEN)<(XLEN-2) && grad[k+XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[i]-1] += 
                    (dzweight + drweight) / sqrt(8.0);
            }
            ++j;
        }
    }
    /*    for(i = nchans;i < toteqn; ++i){
        for(k=0;k < lda;++k){
            printf("%g ",imat[i*lda+k]);
        }
        printf("\n");
    }*/


    eqn = 0;
    /*u = calloc(sizeof(PREC)*2,eqn * lda);*/

    /*uz = calloc(sizeof(PREC)*2,eqn);*/
    nma = 0;


    /*printf("constraint matrix is %d x %d\n",numunk,nma);*/

    sol = calloc(sizeof(PREC) * 2,lda);
    tz = calloc(sizeof(PREC) * 2,toteqn);
    for(i = 0; i < nchans; ++i)tz[i] = inproj[chans[i]];

    for(i = 0; i < toteqn ; ++i){
        for(j = 0; j < lda; ++j){
            timat[j * toteqn + i] = imat[i *lda + j];
        }
    }

    lwork = 10*(lda + nma + toteqn);
    work = calloc(sizeof(PREC) * 2,lwork);
    l_m = toteqn;
    l_n = lda;
    l_p = nma;
    l_lda = toteqn;
    l_ldb = eqn;

    mg = lda;
    me = nma;
    ma = toteqn;
    g = calloc(sizeof(float),mg * mg);
    h = calloc(sizeof(float),mg);
    for(j = 0; j < mg; ++j){
        g[j * mg + j] = 1.0;
    }

    mdw = me + ma + mg;
    w = calloc(sizeof(float),mdw * (lda+1));

    if(me > 0){
        for(i = 0; i < me; ++i){
            for(j = 0; j < lda; ++j){
                w[i + j * mdw] = u[i + j * eqn];
            }
        }
        for(i = 0; i < me; ++i){
            w[i + lda * mdw] = uz[i];
        }
    }

    if(ma > 0){
        for(i = 0; i < toteqn; ++i){
            for(j = 0; j < lda; ++j){
                w[i + me + j * mdw] = imat[j + i * lda];
            }
        }
        for(i = 0;i < toteqn; ++i){
            w[i + me + lda * mdw] = tz[i];
        }
    }

    if(mg > 0){
        for(i = 0;i < mg; ++i){
            for(j = 0; j < lda; ++j){
                w[i + me + ma + j * mdw] = g[i * lda + j];
            }
        }
        for(i = 0;i < mg; ++i){
            w[i + me + ma + lda * mdw] = h[i];
        }
    }




    ip = calloc(sizeof(int),IP_LEN);
    ws = calloc(sizeof(float),WS_LEN);
    ip[0] = WS_LEN;
    ip[1] = IP_LEN;
    prgopt[0] = 4;
    prgopt[1] = 4;
    prgopt[2] = 1.0e-4;
    prgopt[3] = 7;
    prgopt[4] = 5;
    prgopt[5] = 1.0e-4;
    prgopt[6] = 9;
    prgopt[7] = 1;
    prgopt[8] = 1;
    LSEI(w,&mdw,&me,&ma,&mg,&lda,prgopt,sol,&rnorme,&rnorml,&mode,ws,ip);
    printf("lsei mode = %d\n",mode);
    memset(fitimage,0,XLEN*YLEN*sizeof(float));
    for(j = 0,k = 0;k < na;++k){
        if(grad[k]){
            fitimage[k] = sol[j];
            ++j;
        }
    }








    /*for(i = 0 ; i < lda; ++i)printf("sol %d %g\n",i,sol[i]);*/

    /*write_imagedat("core.sdt",fitimage);*/

    bolomclip(shot,fitimage,chans,nchans);
    bolomproj(shot,fitimage,fitproj);
    *fret = 0.0;
    for(i = 0; i < nchans; ++i){
        *fret += pow((inproj[chans[i]] - 
            fitproj[chans[i]]),2.0) / pow(sigma[chans[i]],2.0);
    }
    if(nchans > lda)
        *fret = *fret / (float)(nchans - lda);
    /*fprintf(stderr,"core chi = %g\n",*fret);*/
    free(cmat);
    free(imat);
    free(timat);
    free(amat);
    free(sol);
    free(tz);
    free(work);
    free(g);
    free(h);
    free(w);
    free(ws);

}
/*matherr(x)
register struct exception *x;
{
    if(strcmp(x->name,"sqrt") == 0){
        if(x->arg1 < 0.0)x->retval = sqrt(-x->arg1);
        else x->retval = 0;
        return(1);
    } else if(strcmp(x->name,"log") == 0){
        if(x->arg1 < 0.0)x->retval = -MAXFLOAT;
        else if(x->arg1 == 0.0)x->retval = -MAXFLOAT;
        else x->retval = -MAXFLOAT;
        return(1);
    } else if(strcmp(x->name,"pow") == 0){
        if(x->arg1 < 0.0 && (floor(x->arg2) - x->arg2) != 0.0 )
            x->retval = pow(x->arg1,floor(x->arg2));
        else if(x->arg1 == 0.0)x->retval = 0;
        else x->retval = -MAXFLOAT;
        return(1);
    }
} */